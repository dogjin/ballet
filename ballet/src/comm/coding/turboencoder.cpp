#include "turboencoder.h"
#include "turboencoder_p.h"
#include "convolutionalcoding.h"
#include <splib/matfunc.h>

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
        rintrlvrIndices.set_size(InterleaverIndices.length());
        for (size_t ii=0; ii<InterleaverIndices.length(); ii++)
            rintrlvrIndices(InterleaverIndices(ii)) = ii;

        // set property_locked flag
        property_locked = true;
        
    }

    void TurboEncoderPrivate::encode(const splib::ivec &x, splib::ivec *y)
    {

        // pointer to trellis members
        const int * NEXT = TrellisStructure.nextStates.begin();
        const int * OUT = TrellisStructure.codeOutputs.begin();

        // number of trellis states
        int numStates = TrellisStructure.numStates;

        // number of bits in coded message
        int cwlen = x.length();

        // calculate number of uncoded bits
        int L = (2*N-1)*cwlen + 2*numTails;

        // 1:st constituent encoder
        splib::ivec y1(N*cwlen+numTails);
        {

            const int * X = x.begin();

            // encode
            int s = 0;
            int * Y1 = y1.begin();
            for (size_t i=0; i<cwlen; i++)
            {

                // encoder input
                int inpt = X[i];

                // encoder output codeword (given input and state)
                int otpt = OUT[inpt*numStates+s];

                // convert output codeword to bit representation
                for (size_t n=0; n<N; n++)
                {
                    Y1[N-n-1] = otpt&1;
                    otpt = otpt>>1;
                }

                // get next encoder state
                s = NEXT[inpt*numStates+s];

                // update pointer location
                Y1 = Y1 + N;

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
                    if (NEXT[j*numStates+s] == nxt)
                    {
                        inpt = j;
                        break;
                    }
                }

                // find termination bit output (given state and termination bit input)
                int otpt = OUT[inpt*numStates+s];

                // convert output codeword to bit representation
                for (size_t n=0; n<N; n++)
                {
                    Y1[N-n-1] = otpt&1;
                    otpt = otpt>>1;
                }

                // get next encoder state
                s = nxt;

                // update pointer location
                Y1 = Y1 + N;

            }

        }

        // 2:nd constituent encoder
        splib::ivec y2(N*cwlen+numTails);
        {

            // permute input bits
            const splib::ivec xInv = x(InterleaverIndices);
            const int * X = xInv.begin();

            // encode
            int s = 0;
            int * Y2 = y2.begin();
            for (size_t i=0; i<cwlen; i++)
            {

                // encoder input
                int inpt = X[i];

                // encoder output codeword (given input and state)
                int otpt = OUT[inpt*numStates+s];

                // convert output codeword to bit representation
                for (size_t n=0; n<N; n++)
                {
                    Y2[N-n-1] = otpt&1;
                    otpt = otpt>>1;
                }

                // get next encoder state
                s = NEXT[inpt*numStates+s];

                // update pointer location
                Y2 = Y2 + N;

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
                    if (NEXT[j*numStates+s] == nxt)
                    {
                        inpt = j;
                        break;
                    }
                }

                // find termination bit output (given state and termination bit input)
                int otpt = OUT[inpt*numStates+s];

                // convert output codeword to bit representation
                for (size_t n=0; n<N; n++)
                {
                    Y2[N-n-1] = otpt&1;
                    otpt = otpt>>1;
                }

                // get next encoder state
                s = nxt;

                // update pointer location
                Y2 = Y2 + N;

            }

        }

        // mux turbo encoded data
        y->set_size(L);
        {

            // declare pointers
            int * Y = y->begin();
            int * Y1 = y1.begin();
            int * Y2 = y2.begin();

            for (size_t i=0; i<cwlen; i++)
            {

                // copy data to return array
                // (x,y1,y2)
                Y[0] = Y1[0];
                for (size_t n=0; n<N-1; n++)
                    Y[1+n] = Y1[1+n];
                for (size_t n=0; n<N-1; n++)
                    Y[N+n] = Y2[1+n];

                // update pointers
                Y = Y + 2*N - 1;
                Y1 = Y1 + N;
                Y2 = Y2 + N;

            }

            // [termination bits for 1:st constituent encoder]
            for (size_t i=0; i<K; i++)
            {
                for (size_t n=0; n<N; n++)
                    Y[n] = Y1[n];
                Y = Y + N;
                Y1 = Y1 + N;
            }

            // [termination bits for 2:nd constituent encoder]
            for (size_t i=0; i<K; i++)
            {
                for (size_t n=0; n<N; n++)
                    Y[n] = Y2[n];
                Y = Y + N;
                Y2 = Y2 + N;
            }

        }

    }
    
    TurboEncoder::TurboEncoder()
        : balletObject(*new TurboEncoderPrivate)
    {

        // default trellis
        splib::ivec ConstraintLength = "4";
        splib::imat CodeGenerator = "13 15";
        splib::ivec FeedbackConnection = "13";
        TrellisStructure = ballet::poly2trellis(
            ConstraintLength,CodeGenerator,FeedbackConnection);

        // default InterleaverIndices
        InterleaverIndices = to_ivec(splib::linspace(63,1,64));

    }

    splib::ivec TurboEncoder::encode(const splib::ivec &x)
    {

        BALLET_D(TurboEncoder);

        // lock Encoder object
        if (!isLocked()) { d->lock(); }

        splib::ivec y;
        d->encode(x,&y);

        return y;

    }

};
