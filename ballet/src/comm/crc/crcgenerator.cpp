#include "crcbase_p.h"
#include "crcgenerator.h"

using namespace arma;

namespace ballet
{

    class CRCGeneratorPrivate : public CRCBasePrivate
    {
    public:
        imat generate(const imat &x);
    };

    imat CRCGeneratorPrivate::generate(const imat &x)
    {
        
        // call the baseclass generate function
        imat encData = baseGenerate(x);

        return encData;

    }

    CRCGenerator::CRCGenerator()
        : CRCBase(*new CRCGeneratorPrivate)
    {}

    CRCGenerator::CRCGenerator(unsigned int poly)
        : CRCBase(*new CRCGeneratorPrivate,poly)
    {}

    imat CRCGenerator::generate(const imat &x)
    {
        
        BALLET_D(CRCGenerator);
        if (!isLocked()) { d->lock(); }
        return d->generate(x);

    }
        

};
