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
        splib::Trellis TrellisStructure;
        splib::ivec InterleaverIndices;

    public:
        int N;
        int K;
        int numTails;

    public:
        splib::ivec rintrlvrIndices;

    public:
        void encode(const splib::ivec &x, splib::ivec *y);
        void lock();

    };

};

#endif
