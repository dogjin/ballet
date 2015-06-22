#ifndef __ballet_viterbi_p_h
#define __ballet_viterbi_p_h

// ballet includes
#include "viterbi.h"
#include "private/object_p.h"

namespace ballet
{

    struct VDExtra
    {


    };

    class ViterbiDecoderPrivate : public balletObjectPrivate
    {

        BALLET_DECLARE_PUBLIC(ViterbiDecoder)

    public:
        ViterbiDecoderPrivate();
        ~ViterbiDecoderPrivate();

    public:
        Trellis TrellisStructure;
        std::string InputFormat;
        size_t SoftInputWordLength;
        std::string OutputFormat;

    public:
        int numTails;

    public:
        arma::imat decode(const arma::mat &x);
        void lock();

    public:
        size_t k;           // number of bits per input symbol
        size_t n;           // number of bits per output symbol
        
        double soft_zero;   // most confident logical zero
        double soft_one;    // most confident logical one

        arma::mat binout;   // output codeword bits

    protected:
        arma::mat branch_metric_unq(const arma::mat &x);
        arma::mat branch_metric_har(const arma::mat &x);
        arma::mat branch_metric_sof(const arma::mat &x);

    };

};

#endif
