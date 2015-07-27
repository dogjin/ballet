// ballet includes
#include "turboencoder.h"
#include "turboencoder_p.h"
#include "convolutionalcoding.h"

using namespace arma;

namespace ballet
{

    TurboEncoderPrivate::TurboEncoderPrivate()
        : balletObjectPrivate()
    {
    }

    TurboEncoderPrivate::~TurboEncoderPrivate()
    {
    }

    void TurboEncoderPrivate::lock()
    {

        BALLET_Q(TurboEncoder);

        // lock TrellisStructure property
        TrellisStructure = q->TrellisStructure;

        // lock InputFormat property
        InterleaverIndices = q->InterleaverIndices;

        // calculate encoder memory
        double Kdbl = log((double)TrellisStructure.numStates)/log(2.0);
        K = static_cast<int>(Kdbl+0.5);
	
        // calculate number of bits per input symbol
        double Ndbl = log((double)TrellisStructure.numOutputSymbols)/log(2.0);
        N = static_cast<int>(Ndbl+0.5);

        // calculate number of tail bits
        numTails = K * N;

        // map reverse interleaver
        rintrlvrIndices.set_size(
            InterleaverIndices.n_rows,
            InterleaverIndices.n_cols);
        for (size_t ii=0; ii<InterleaverIndices.n_elem; ii++)
            rintrlvrIndices(InterleaverIndices(ii)) = ii;

        // set property_locked flag
        property_locked = true;
        
    }

    void TurboEncoderPrivate::encode(const imat &x, imat *y)
    {

        // number of trellis states
        int numStates = TrellisStructure.numStates;

        // number of bits in coded message
        int cwlen = x.n_elem;

        // calculate number of uncoded bits
        int L = (2*N-1)*cwlen + 2*numTails;

        // 1:st constituent encoder
        imat y1(N*cwlen+numTails,1);
        {

            const int * X = x.begin();

            // encode
            int s = 0;
            for (size_t i=0; i<cwlen; i++)
            {

                // encoder input
                int inpt = x(i);

                // encoder output codeword (given input and state)
                int otpt = TrellisStructure.outputs(s,inpt);

                // convert output codeword to bit representation
                for (size_t n=0; n<N; n++)
                {
                    y1(N*(i+1)-n-1) = otpt&1;
                    otpt = otpt>>1;
                }

                // get next encoder state
                s = TrellisStructure.nextStates(s,inpt);

            }

            // [trellis termination]
            for (size_t i=0; i<K; i++)
            {

                // peek next state
                int nxt = s >> 1;

                // find termination bit
                int inpt = -1;
                for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                {
                    if (TrellisStructure.nextStates(s,j) == nxt)
                    {
                        inpt = j;
                        break;
                    }
                }

                // find termination bit output (given state and termination bit input)
                int otpt = TrellisStructure.outputs(s,inpt);

                // convert output codeword to bit representation
                for (size_t n=0; n<N; n++)
                {
                    y1(N*(cwlen+i+1)-n-1) = otpt&1;
                    otpt = otpt>>1;
                }

                // get next encoder state
                s = nxt;

            }

        }

        // 2:nd constituent encoder
        imat y2(N*cwlen+numTails,1);
        {

            // permute input bits
            const imat xInv = x.elem(InterleaverIndices);

            // initial state
            int s = 0;

            // encode
            for (size_t i=0; i<cwlen; i++)
            {

                // encoder input
                int inpt = xInv(i);

                // encoder output codeword (given input and state)
                int otpt = TrellisStructure.outputs(s,inpt);

                // convert output codeword to bit representation
                for (size_t n=0; n<N; n++)
                {
                    y2(N*(i+1)-n-1) = otpt&1;
                    otpt = otpt>>1;
                }

                // get next encoder state
                s = TrellisStructure.nextStates(s,inpt);

            }

            // [trellis termination]
            for (size_t i=0; i<K; i++)
            {

                // peek next state
                int nxt = s >> 1;

                // find termination bit
                int inpt = -1;
                for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                {
                    if (TrellisStructure.nextStates(s,j) == nxt)
                    {
                        inpt = j;
                        break;
                    }
                }

                // find termination bit output (given state and termination bit input)
                int otpt = TrellisStructure.outputs(s,inpt);

                // convert output codeword to bit representation
                for (size_t n=0; n<N; n++)
                {
                    y2(N*(cwlen+i+1)-n-1) = otpt&1;
                    otpt = otpt>>1;
                }

                // get next encoder state
                s = nxt;

            }

        }

        // mux turbo encoded data
        y->set_size(L,1);
        {

            // convert object pointer to reference object
            imat & tmp001 = *y;

            // declare pointers
            int * Y = y->begin();
            int * Y1 = y1.begin();
            int * Y2 = y2.begin();

            for (size_t i=0; i<cwlen; i++)
            {

                // turbo encoder interleaver
                tmp001((2*N-1)*i) = y1(N*i);
                for (size_t n=0; n<N-1; n++)
                    tmp001((2*N-1)*i+n+1) = y1(N*i+n+1);
                for (size_t n=0; n<N-1; n++)
                    tmp001((2*N-1)*i+n+N) = y2(N*i+n+1);

            }

            // [termination bits for 1:st constituent encoder]
            for (size_t i=0; i<K; i++)
            {
                for (size_t n=0; n<N; n++)
                    tmp001((2*N-1)*cwlen+N*i+n) = y1(N*cwlen+N*i+n);
            }

            // [termination bits for 2:nd constituent encoder]
            for (size_t i=0; i<K; i++)
            {
                for (size_t n=0; n<N; n++)
                    tmp001((2*N-1)*cwlen+N*K+N*i+n) = y2(N*cwlen+N*i+n);
            }

        }

    }
    
    TurboEncoder::TurboEncoder()
        : balletObject(*new TurboEncoderPrivate)
    {

        // default trellis
        imat ConstraintLength = "4";
        imat CodeGenerator = "13 15";
        imat FeedbackConnection = "13";
        TrellisStructure = poly2trellis(
            ConstraintLength,CodeGenerator,FeedbackConnection);

        // default InterleaverIndices
        InterleaverIndices = linspace<umat>(63,0,64);

    }

    imat TurboEncoder::encode(const imat &x)
    {

        BALLET_D(TurboEncoder);

        // lock Encoder object
        if (!isLocked()) { d->lock(); }

        imat y;
        d->encode(x,&y);

        return y;

    }

};
