#ifndef __ballet_crcdetector_h
#define __ballet_crcdetector_h

#include "crcbase.h"

namespace ballet
{

    class CRCDetectorPrivate;

    class CRCDetector : public CRCBase
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(CRCDetector)
        /// @endcond

    public:
        CRCDetector();
        CRCDetector(unsigned int poly);

    public:
        int detect(const arma::imat &x);

    };

};

#endif
