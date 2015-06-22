#ifndef __ballet_comm_turboencoder_h
#define __ballet_comm_turboencoder_h

#include <vector>
#include <splib/vec.h>
#include <splib/commfunc.h>
#include "ballet/object.h"

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
        splib::Trellis TrellisStructure;
        splib::ivec InterleaverIndices;

    public:
        splib::ivec encode(const splib::ivec &x);

    };

};

#endif
