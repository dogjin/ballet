#ifndef __ballet_comm_de2bi_h
#define __ballet_comm_de2bi_h

#include <armadillo>

namespace ballet
{

    arma::imat de2bi(const int d);
    arma::imat de2bi(const int d, const int n);
    arma::imat de2bi(const arma::imat d);
    arma::imat de2bi(const arma::imat d, const int n);

};

#endif
