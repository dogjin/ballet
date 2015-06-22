#ifndef __ballet_crcbase_p_h
#define __ballet_crcbase_p_h

#include "crcbase.h"
#include "private/object_p.h"

namespace ballet
{

    class CRCBasePrivate : public balletObjectPrivate
    {

        BALLET_DECLARE_PUBLIC(CRCBase)

    public:
        CRCBasePrivate();

    public:
        void lock();

    public:
        arma::imat Polynomial;
        bool ReflectChecksums;

    public:
        arma::imat baseGenerate(const arma::imat &x);

    public:
        unsigned int crc_size;
        unsigned int generator;
        unsigned int crc_mask;
        unsigned int crc_topbit;
        void initObject();

    };

};

#endif
