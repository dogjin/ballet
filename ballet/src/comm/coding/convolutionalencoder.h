#ifndef __ballet_comm_convolutionalencoder_h
#define __ballet_comm_convolutionalencoder_h

// ballet includes
#include "ballet/object.h"
#include "convolutionalcoding.h"

// external library includes
#include <armadillo>

namespace ballet
{

    class ConvolutionalEncoderPrivate;

    class ConvolutionalEncoder : public balletObject
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(ConvolutionalEncoder)
        /// @endcond

    public:
        ConvolutionalEncoder();

    public:
        Trellis TrellisStructure;

    public:
        arma::imat encode(const arma::imat &x);

    };

};


#endif
