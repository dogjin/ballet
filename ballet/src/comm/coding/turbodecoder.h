#ifndef __ballet_comm_turbodecoder_h
#define __ballet_comm_turbodecoder_h

#include <vector>
#include <splib/vec.h>
#include <splib/commfunc.h>
#include "ballet/object.h"

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
        splib::Trellis TrellisStructure;
        splib::ivec InterleaverIndices;
        std::string Algorithm;
        int NumScalingBits;
        int NumIterations;

    public:
        splib::fvec decode(const splib::fvec &x);

    };

};

#endif
