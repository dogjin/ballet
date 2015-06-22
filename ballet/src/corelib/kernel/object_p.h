#ifndef __ballet_object_p_h
#define __ballet_object_p_h

#include "object.h"

namespace ballet
{

    class balletObjectPrivate : public balletObjectData
    {

        BALLET_DECLARE_PUBLIC(balletObject)

    public:
        balletObjectPrivate();
        virtual ~balletObjectPrivate();

    public:
        bool isLocked();
        void lock();
        void release();

    };

};

#endif
