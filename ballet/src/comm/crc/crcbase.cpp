#include "crcbase_p.h"
#include <ballet/comm/de2bi.h>

using namespace arma;

namespace ballet
{

    imat CRCBase::convertHexToBinary(unsigned int hexval)
    {
        
        // convert decimal base to hexidecimal base
        double hexBase = std::log((double)hexval) / std::log(16);

        // determine number of hexidecimal places
        size_t lenHexVal = static_cast<size_t>(std::ceil(hexBase));
        if (hexval >= (unsigned int)(1<<(4*lenHexVal))) { lenHexVal++; }

        // preallocate binary vector
        imat binaryVector(4*lenHexVal,1);

        // convert to binary by one hex symbol at a time
        for (size_t idx=0; idx<lenHexVal; idx++)
        {

            // get hex symbol at idx
            unsigned int hexSymbol = (hexval>>(4*(lenHexVal-idx-1))) & 0xf;

            // convert hext symbol at idx to binary
            for (size_t b=0; b<4; b++)
            {
                binaryVector(4*idx+3-b) = hexSymbol & 1;
                hexSymbol = hexSymbol >> 1;
            }

        }

        return binaryVector;

    }

    CRCBasePrivate::CRCBasePrivate()
        : balletObjectPrivate()
    {}

    void CRCBasePrivate::lock()
    {

        BALLET_Q(CRCBase);

        // lock Polynomial property
        Polynomial = q->Polynomial;

        // lock ReflectChecksums property
        ReflectChecksums = q->ReflectChecksums;

        // initialize crc object
        initObject();

        // set property_locked flag
        property_locked = true;

    }

    /*
      Initialize CRCBase object
    */
    void CRCBasePrivate::initObject()
    {

        // get the length of the polynomial
        imat poly = Polynomial(span(1,Polynomial.n_rows-1),0);
        crc_size = poly.n_elem;

        // translate the generator polynomail to hexidecimal
        generator = 0x0;
        for (size_t ii=0; ii<crc_size;  ii++)
        {
            generator = generator | (unsigned int)(poly(ii)<<(crc_size-ii-1));
        }

        // checksum mask
        crc_mask = ((unsigned int)(1<<(crc_size-1))-1)<<1|1;

        // checksum topmost bit
        crc_topbit = 1<<(crc_size-1);

    }

    imat CRCBasePrivate::baseGenerate(const imat &x)
    {

        // codeword size
        unsigned int cwsize = x.n_elem;

        // get the length of the polynomial
        imat poly = Polynomial(span(1,Polynomial.n_rows-1),0);
        int polyLen = poly.n_elem;

        // Now proceed as if it's a regular CRC computation
        // Polynomial division with deferred message XORing
        unsigned int remainder = 0x0; 
        for (int i=0; i<cwsize; i++)
        {

            remainder = remainder ^ (x(i) << (crc_size-1));
            if ( remainder & crc_topbit )
                remainder = (remainder<<1) ^ generator;
            else
                remainder = remainder << 1;

        }

        // apply crc mask
        remainder = remainder & crc_mask;

        // convert remainder from decimal
        // to binary representation
        imat reg = de2bi(remainder,polyLen);

        // reflect the checksum if appropriate
        if (ReflectChecksums)
        {
            reg = flipud(reg);
        }

        // add checksum
        imat encData(cwsize+polyLen,1);
	encData(span(0,cwsize-1),0) = x;
	encData(span(cwsize,cwsize+polyLen-1),0) = reg;

        return encData;

    }

    CRCBase::CRCBase()
        : balletObject(*new CRCBasePrivate)
    {

        // initialize default object
        reset();

    }

    CRCBase::CRCBase(unsigned int poly)
        : balletObject(*new CRCBasePrivate)
    {

         // initialize default object
        reset();

        // set polynomial
        imat polynomial = convertHexToBinary(poly);
        Polynomial.resize(1+polynomial.n_rows,1);
        Polynomial(0) = 1;
        Polynomial(span(1,polynomial.n_rows),0) = polynomial;

    }

    CRCBase::CRCBase(CRCBasePrivate &dd)
        : balletObject(dd)
    {

         // initialize default object
        reset();

    }

     CRCBase::CRCBase(CRCBasePrivate &dd, unsigned int poly)
        : balletObject(dd)
    {

         // initialize default object
        reset();

        // set polynomial
        imat polynomial = convertHexToBinary(poly);
        Polynomial.resize(1+polynomial.n_rows,1);
        Polynomial(0) = 1;
        Polynomial(span(1,polynomial.n_rows),0) = polynomial;

    }

    CRCBase::~CRCBase() {}

    void CRCBase::reset()
    {

        // default polynomial = 0x1021
        Polynomial = "1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1";

        // default ReflectChecksums = false
        ReflectChecksums = false;

        // unlock object [if applicable]
        if (isLocked()) { release(); }

    }

    imat CRCBase::baseGenerate(const imat &x)
    {

        BALLET_D(CRCBase);

        // lock crc object
        if (!isLocked()) { d->lock(); }

        return d->baseGenerate(x);

    }

};
