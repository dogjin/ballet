#ifndef __ballet_turbodecoder_p_h
#define __ballet_turbodecoder_p_h

#include "turbodecoder.h"
#include "appdecoder.h"
#include "private/object_p.h"

namespace ballet
{

    class TurboDecoderPrivate : public balletObjectPrivate
    {

        BALLET_DECLARE_PUBLIC(TurboDecoder)

    public:
        TurboDecoderPrivate();
        ~TurboDecoderPrivate();

    public:
        splib::Trellis TrellisStructure;
        splib::ivec InterleaverIndices;
        std::string Algorithm;
        int NumScalingBits;
        int NumIterations;

    public:
        ballet::APPDecoder decObj;

    public:
        int N;
        int K;
        int numTails;

    public:
        splib::ivec rintrlvrIndices;

    public:
        void decode(const splib::fvec &x, splib::fvec *y);
        void lock();

    };

};

#endif
