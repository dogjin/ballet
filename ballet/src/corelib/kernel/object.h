#ifndef __ballet_object_h
#define __ballet_object_h

#include <ballet/global.h>

namespace ballet
{

    class balletObject;
    class balletObjectPrivate;

    class balletObjectData
    {
    public:
        virtual ~balletObjectData() = 0;
        balletObject *q_ptr;
    public:
        bool property_locked;
    };

    class balletObject
    {

        BALLET_DECLARE_PRIVATE(balletObject)

    public:
        balletObject();
        ~balletObject();

    public:
        bool isLocked();
        void release();

    protected:
        balletObject(balletObjectPrivate &dd);

    protected:
        balletObjectData *d_ptr;
    };

};

#endif
