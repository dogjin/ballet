#ifndef __ballet_comm_turbodecoder_h
#define __ballet_comm_turbodecoder_h

#include <armadillo>
#include "ballet/object.h"
#include "convolutionalcoding.h"

namespace ballet
{

    class TurboDecoderPrivate;

    class TurboDecoder : public balletObject
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(TurboDecoder)
        /// @endcond

    public:
        TurboDecoder();

    public:
        Trellis TrellisStructure;
        arma::uvec InterleaverIndices;
        std::string Algorithm;
        int NumScalingBits;
        int NumIterations;

    public:
        arma::mat decode(const arma::mat &x);

    };

};

#endif
