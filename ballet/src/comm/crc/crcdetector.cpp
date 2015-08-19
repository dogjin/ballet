#include "crcbase_p.h"
#include "crcdetector.h"

using namespace arma;

namespace ballet
{

    class CRCDetectorPrivate : public CRCBasePrivate
    {
    public:
        int detect(const imat &x);
    };

    int CRCDetectorPrivate::detect(const imat &x)
    {

        // get the length of the polynomial
        imat poly = Polynomial(0,span(1,Polynomial.n_cols-1));
        int polyLen = poly.n_elem;

        // get the length of the input codeword
        unsigned int cwsize = x.n_elem - polyLen;

        // call the baseclass generate function on the
        // original data, then compare its output with
        // the checksum
        imat outData = x(span(0,cwsize-1),0);
        imat encData = baseGenerate(outData);

        // determine if an error has occured
        imat inChecksum = x(span(cwsize,cwsize+polyLen-1),0);
        imat outChecksum = encData(span(cwsize,cwsize+polyLen-1),0);
        bool error = any(vectorise(inChecksum != outChecksum));

        return static_cast<int>(error);

    }

    CRCDetector::CRCDetector()
        : CRCBase(*new CRCDetectorPrivate)
    {}

    CRCDetector::CRCDetector(unsigned int poly)
        : CRCBase(*new CRCDetectorPrivate,poly)
    {}

    int CRCDetector::detect(const imat &x)
    {
        
        BALLET_D(CRCDetector);
        if (!isLocked()) { d->lock(); }
        return d->detect(x);

    }
        

};
