#ifndef __ballet_crcgenerator_h
#define __ballet_crcgenerator_h

#include "crcbase.h"

namespace ballet
{

    class CRCGeneratorPrivate;

    class CRCGenerator : public CRCBase
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(CRCGenerator)
        /// @endcond

    public:
        CRCGenerator();
        CRCGenerator(unsigned int poly);

    public:
        arma::imat generate(const arma::imat &x);

    };

};

#endif
