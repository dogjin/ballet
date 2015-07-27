#ifndef __ballet_comm_turboencoder_h
#define __ballet_comm_turboencoder_h

// ballet includes
#include "ballet/object.h"
#include "convolutionalcoding.h"

// external includes
#include <armadillo>

namespace ballet
{

    class TurboEncoderPrivate;

    class TurboEncoder : public balletObject
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(TurboEncoder)
        /// @endcond

    public:
        TurboEncoder();

    public:
        Trellis TrellisStructure;
        arma::umat InterleaverIndices;

    public:
        arma::imat encode(const arma::imat &x);

    };

};

#endif
