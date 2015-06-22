#include "object.h"
#include "object_p.h"

namespace ballet
{

    balletObjectData::~balletObjectData() {}

    balletObjectPrivate::balletObjectPrivate()
    {

        // balletObjectData initialization
        q_ptr = 0;
        property_locked = false;

    }

    bool balletObjectPrivate::isLocked()
    {
        return property_locked;
    }

    void balletObjectPrivate::lock()
    {
        property_locked = true;
    }

    void balletObjectPrivate::release()
    {
        property_locked = false;
    }

    balletObjectPrivate::~balletObjectPrivate() {}

    balletObject::balletObject()
        : d_ptr(new balletObjectPrivate)
    {
        BALLET_D(balletObject);
        d_ptr->q_ptr = this;
    }

    /*!
        \internal
    */
    balletObject::balletObject(balletObjectPrivate &dd)
        : d_ptr(&dd)
    {
        BALLET_D(balletObject);
        d_ptr->q_ptr = this;
    }

    balletObject::~balletObject() { delete d_ptr; }

    bool balletObject::isLocked()
    {
        BALLET_D(balletObject);
        return d->isLocked();
    }

    void balletObject::release()
    {
        BALLET_D(balletObject);
        d->release();
    }

};
