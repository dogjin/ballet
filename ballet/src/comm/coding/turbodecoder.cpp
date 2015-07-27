#include "turbodecoder.h"
#include "turbodecoder_p.h"
#include <ballet/logical.h>

using namespace arma;

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
        rintrlvrIndices.copy_size(InterleaverIndices);
        for (size_t ii=0; ii<InterleaverIndices.n_elem; ++ii)
        {
            rintrlvrIndices(InterleaverIndices(ii)) = ii;
        }

        // setup APPDecoder object
        decObj.TrellisStructure = TrellisStructure;
        decObj.Algorithm = Algorithm;
        decObj.NumScalingBits = NumScalingBits;

        // set property_locked flag
        property_locked = true;
        
    }

    void TurboDecoderPrivate::decode(const mat &x, mat *_y)
    {

        // number of bits in coded message
        int cwlen = x.n_elem;

        // calculate number of uncoded bits
        int L = (cwlen-2*numTails) / (2*N-1);

        // get systematic bits
        vec xRx(L+K);
        vec xIn(L+K);
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
        vec lc1(N*(L+K));
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
        vec lc2(N*(L+K));
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
        vec Le12 = zeros<vec>(L+K);
        vec Le21 = zeros<vec>(L+K);
        vec Le;

        // declare decoder result vector
        ivec dec;

        // declare and initialize comparison vector
        ivec prior = zeros<ivec>(L);

        // declare and initialize zero vector
        vec myzeros = zeros<vec>(3);

        int iteration = 0;
        while (iteration < NumIterations)
        {

            // 1:st constituent decoder
            vec La1 = Le21.elem(rintrlvrIndices);
            La1.insert_rows(La1.n_elem,3);
            vec LU1 = decObj.decode(La1,lc1);
            Le12 = LU1;

            // 2:nd constituent decoder
            vec La2 = Le12.elem(InterleaverIndices);
            La2.insert_rows(La2.n_elem,3);
            vec LU2 = decObj.decode(La2,lc2);
            Le21 = LU2;

            // convert to bit-decisions
            Le = Le12(span(0,L-1)) + Le21.elem(rintrlvrIndices);
            dec = getHardDecisions(Le);

            if (ballet::all(dec == prior)) { break; }
            prior = dec;

            iteration++;

        }

        // resize return vector
        mat& y = *_y;
        y = zeros<mat>(L,1);
    
        // copy result to return vector
        for (size_t i=0; i<L; i++)
            y(i) = static_cast<double>(dec(i));

    }
    
    TurboDecoder::TurboDecoder()
        : balletObject(*new TurboDecoderPrivate)
    {

        // default trellis
        imat ConstraintLength = "4";
        imat CodeGenerator = "13 15";
        imat FeedbackConnection = "13";
        TrellisStructure = ballet::poly2trellis(
            ConstraintLength,CodeGenerator,FeedbackConnection);

        // default InterleaverIndices
        InterleaverIndices = linspace<umat>(63,0,64);

        // defalut Algorithm
        Algorithm = "True APP";

        // default NumScalingBits
        NumScalingBits = 3;
        
        // default NumIterations
        NumIterations = 6;

    }

    mat TurboDecoder::decode(const mat &x)
    {

        BALLET_D(TurboDecoder);

        // lock decoder object
        if (!isLocked()) { d->lock(); }

        mat y;
        d->decode(x,&y);

        return y;

    }

};
