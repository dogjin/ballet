// ballet includes
#include "convolutionalencoder.h"
#include "convolutionalencoder_p.h"

// system includes
#include <limits>

using namespace arma;

namespace ballet
{

    ConvolutionalEncoderPrivate::ConvolutionalEncoderPrivate()
        : balletObjectPrivate()
    {
    }

    ConvolutionalEncoderPrivate::~ConvolutionalEncoderPrivate()
    {
    }

    void ConvolutionalEncoderPrivate::lock()
    {
        BALLET_Q(ConvolutionalEncoder);

        // lock TrellisStructure property
        TrellisStructure = q->TrellisStructure;

        // number of bits per symbol
        double kdbl = std::log((double)TrellisStructure.numInputSymbols) / std::log(2.0);
        k = static_cast<size_t>(kdbl);
        double ndbl = std::log((double)TrellisStructure.numOutputSymbols) / std::log(2.0);
        n = static_cast<size_t>(ndbl);

        // number of tail bits
        double numTails_ = log((double)TrellisStructure.numStates)/log(2.0);
        numTails = static_cast<int>(numTails_+0.5);

        // set property_locked flag
        property_locked = true;
        
    }

    imat ConvolutionalEncoderPrivate::encode(const imat &x)
    {

        // large floating point number
        double inf = std::numeric_limits<double>::max();

        // number of trellis states
        size_t numStates = TrellisStructure.numStates;

        // uncoded frame length
        size_t frmLen = x.n_elem;

        // number number of trellis sections
        size_t L = frmLen/k;
       
        // coded frame length
        size_t cdeLen = (L+numTails) * n;

        // initialize output array
        imat y(cdeLen,1,fill::zeros);

        // set initial state 
        int s = 0;

        for (size_t ii=1; ii<=L; ii++)
        {

            // compute uncoded input word from input bits
            int input = 0;
            for (size_t kk=0; kk<k; ++kk)
            {
                input = input | (x(ii*k-kk-1)<<kk);
            }

            // compute encoder output codeword
            int output = TrellisStructure.outputs(s,input);

            // convert output codeword to binary
            for (size_t nn=0; nn<n; ++nn)
            {
                if (output == 0) { break; }
                y(ii*n-nn-1) = output&1;
                output = output>>1;
            }

            // compute next trellis state
            s = TrellisStructure.nextStates(s,input);

        }

        for (size_t ii=L+1; ii<=L+numTails; ii++)
        {

            // compute encoder output codeword
            int output = TrellisStructure.outputs(s,0);

            // convert output codeword to binary
            for (size_t nn=0; nn<n; nn++)
            {
                if (output == 0) { break; }
                y(ii*n-nn-1) = output&1;
                output = output>>1;
            }

            // compute next trellis state
            s = TrellisStructure.nextStates(s,0);

        }

        return y;

    }
    
    ConvolutionalEncoder::ConvolutionalEncoder()
        : balletObject(*new ConvolutionalEncoderPrivate)
    {

        // default trellis
        imat constraints = "7";
        imat codegen = "177 133;";
        TrellisStructure = poly2trellis(constraints,codegen);

    }

    imat ConvolutionalEncoder::encode(const imat &x)
    {

        BALLET_D(ConvolutionalEncoder);

        // lock encoder object
        if (!isLocked()) { d->lock(); }

        return d->encode(x);

    }

};
