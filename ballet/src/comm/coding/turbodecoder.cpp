#include "turbodecoder.h"
#include "turbodecoder_p.h"
#include "convolutionalcoding.h"
#include <ballet/logical.h>
#include <splib/matfunc.h>

namespace ballet
{

    TurboDecoderPrivate::TurboDecoderPrivate()
        : balletObjectPrivate()
    {
    }

    TurboDecoderPrivate::~TurboDecoderPrivate()
    {
    }

    void TurboDecoderPrivate::lock()
    {

        BALLET_Q(TurboDecoder);

        // lock TrellisStructure property
        TrellisStructure = q->TrellisStructure;

        // lock InputFormat property
        InterleaverIndices = q->InterleaverIndices;
	
        // lock Algorithm property
        Algorithm = q->Algorithm;

        // lock NumScalingBits property
        NumScalingBits = q->NumScalingBits;

        // lock NumIterations property
        NumIterations = q->NumIterations;

        // calculate code rate
        double Ndbl = log((double)TrellisStructure.numOutputSymbols)/log(2.0);
        N = static_cast<int>(Ndbl+0.5);

        // calculate encoder memory
        double Kdbl = log((double)TrellisStructure.numStates)/log(2.0);
        K = static_cast<int>(Kdbl+0.5);

        // calculate number of tail bits
        numTails = K * N;

        // map reverse interleaver
        rintrlvrIndices.set_size(InterleaverIndices.length());
        for (size_t ii=0; ii<InterleaverIndices.length(); ii++)
            rintrlvrIndices(InterleaverIndices(ii)) = ii;

        // setup APPDecoder object
        decObj.TrellisStructure = TrellisStructure;
        decObj.Algorithm = Algorithm;
        decObj.NumScalingBits = NumScalingBits;

        // set property_locked flag
        property_locked = true;
        
    }

    void TurboDecoderPrivate::decode(const splib::fvec &x, splib::fvec *y)
    {

        // number of bits in coded message
        int cwlen = x.length();

        // calculate number of uncoded bits
        int L = (cwlen-2*numTails) / (2*N-1);

        // get systematic bits
        splib::fvec xRx(L+K);
        splib::fvec xIn(L+K);
        for (size_t i=0; i<L; i++)
        {
            xRx(i) = x((2*N-1)*i);
            xIn(i) = 0.0;
        }
        for (size_t i=0; i<K; i++)
        {
            xRx(L+i) = x((2*N-1)*L+N*i);
            xIn(L+i) = x((2*N-1)*L+N*i+numTails);
        }

        // get coded likelihood vector
        // from 1:st constituent encoder
        splib::fvec lc1(N*(L+K));
        for (size_t i=0; i<L; i++)
        {
            lc1(N*i) = x((2*N-1)*i);
            for (size_t n=0; n<N-1; n++)
            {
                lc1(N*i+n+1) = x((2*N-1)*i+n+1);
            }
        }
        for (size_t i=0; i<K; i++)
        {
            lc1(N*(L+i)) = x((2*N-1)*L+N*i);
            for (size_t n=0; n<N-1; n++)
            {
                lc1(N*(L+i)+n+1) = x((2*N-1)*L+N*i+n+1);
            }
        }

        // get coded likelihood vector
        // from 2:nd constituent encoder
        splib::fvec lc2(N*(L+K));
        for (size_t i=0; i<L; i++)
        {
            lc2(N*i) = 0.0;
            for (size_t n=0; n<N-1; n++)
            {
                lc2(N*i+n+1) = x((2*N-1)*i+N+n);
            }
        }
        for (size_t i=0; i<K; i++)
        {
            lc2(N*(L+i)) = x((2*N-1)*L+N*i+numTails);
            for (size_t n=0; n<N-1; n++)
            {
                lc2(N*(L+i)+n+1) = x((2*N-1)*L+N*i+n+1+numTails);
            }
        }

        // declare extrinsic likelihood vectors
        splib::fvec Le12(L+K); Le12.zeros();
        splib::fvec Le21(L+K); Le21.zeros();
        splib::fvec Le;

        // declare decoder result vector
        splib::ivec dec;

        // declare and initialize comparison vector
        splib::ivec prior(L); prior.zeros();

        // declare and initialize zero vector
        splib::fvec myzeros(3);
        myzeros.zeros();

        int iteration = 0;
        while (iteration < NumIterations)
        {

            // 1:st constituent decoder
            splib::fvec La1 = splib::concat(Le21(rintrlvrIndices),myzeros);
            splib::fvec LU1 = decObj.decode(La1,lc1);
            Le12 = LU1;

            // 2:nd constituent decoder
            splib::fvec La2 = splib::concat(Le12(InterleaverIndices),myzeros);
            splib::fvec LU2 = decObj.decode(La2,lc2);
            Le21 = LU2;

            // convert to bit-decisions
            Le = Le12(0,L-1) + Le21(rintrlvrIndices);
            dec = getHardDecisions(Le);

            if (ballet::all(dec == prior)) { break; }
            prior = dec;

            iteration++;

        }

        // resize return vector
        y->set_size(L);
    
        // copy result to return vector
        float * Y = y->begin();
        for (size_t i=0; i<L; i++)
            Y[i] = static_cast<float>(dec[i]);

    }
    
    TurboDecoder::TurboDecoder()
        : balletObject(*new TurboDecoderPrivate)
    {

        // default trellis
        splib::ivec ConstraintLength = "4";
        splib::imat CodeGenerator = "13 15";
        splib::ivec FeedbackConnection = "13";
        TrellisStructure = ballet::poly2trellis(
            ConstraintLength,CodeGenerator,FeedbackConnection);

        // default InterleaverIndices
        InterleaverIndices = to_ivec(splib::linspace(63,1,64));

        // defalut Algorithm
        Algorithm = "True APP";

        // default NumScalingBits
        NumScalingBits = 3;
        
        // default NumIterations
        NumIterations = 6;

    }

    splib::fvec TurboDecoder::decode(const splib::fvec &x)
    {

        BALLET_D(TurboDecoder);

        // lock decoder object
        if (!isLocked()) { d->lock(); }

        splib::fvec y;
        d->decode(x,&y);

        return y;

    }

};
