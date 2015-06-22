#ifndef __ballet_crcbase_h
#define __ballet_crcbase_h

#include <armadillo>
#include "ballet/object.h"

namespace ballet
{

    class CRCBasePrivate;

    class CRCBase : public balletObject
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(CRCBase)
        /// @endcond

    public:
        CRCBase();
        CRCBase(unsigned int poly);
        ~CRCBase();

    protected:
        CRCBase(CRCBasePrivate &dd);
        CRCBase(CRCBasePrivate &dd, unsigned int poly);

    public:
        void reset();

    public:
        arma::imat Polynomial;
        bool ReflectChecksums;

    public:
        arma::imat baseGenerate(const arma::imat &x);

    public:
        static arma::imat convertHexToBinary(unsigned int hexval);

    };

};

#endif
