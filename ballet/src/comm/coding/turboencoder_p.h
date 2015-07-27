#ifndef __ballet_turboencoder_p_h
#define __ballet_turboencoder_p_h

#include "turboencoder.h"
#include "private/object_p.h"

namespace ballet
{

    class TurboEncoderPrivate : public balletObjectPrivate
    {

        BALLET_DECLARE_PUBLIC(TurboEncoder)

    public:
        TurboEncoderPrivate();
        ~TurboEncoderPrivate();

    public:
        Trellis TrellisStructure;
        arma::umat InterleaverIndices;

    public:
        int N;
        int K;
        int numTails;

    public:
        arma::umat rintrlvrIndices;

    public:
        void encode(const arma::imat &x, arma::imat *y);
        void lock();

    };

};

#endif
