#ifndef __ballet_de2bi_p_h
#define __ballet_de2bi_p_h

#include <armadillo>

namespace ballet
{

    class de2biPrivate
    {
    public:
        static arma::imat de2bi(const arma::imat &d, const int n);
    };

};

#endif
