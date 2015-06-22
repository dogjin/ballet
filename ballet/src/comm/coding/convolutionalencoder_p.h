#ifndef __ballet_convolutionalencoder_p_h
#define __ballet_convolutionalencoder_p_h

#include "private/object_p.h"

namespace ballet
{

    class ConvolutionalEncoderPrivate : public balletObjectPrivate
    {

        BALLET_DECLARE_PUBLIC(ConvolutionalEncoder)

    public:
        ConvolutionalEncoderPrivate();
        ~ConvolutionalEncoderPrivate();

    public:
        Trellis TrellisStructure;

    public:
        int numTails;

    public:
        arma::imat encode(const arma::imat &x);
        void lock();

    public:
        size_t k;           // number of bits per input symbol
        size_t n;           // number of bits per output symbol

    protected:
        arma::mat branch_metric_unq(const arma::mat &x);
        arma::mat branch_metric_har(const arma::mat &x);
        arma::mat branch_metric_sof(const arma::mat &x);

    };

};

#endif
