#ifndef __ballet_convolutionalcoding_h
#define __ballet_convolutionalcoding_h

#include <armadillo>

namespace ballet
{

    typedef struct
    {
        int         numInputSymbols;
        int         numOutputSymbols;
        int         numStates;
        arma::imat  nextStates;
        arma::imat  outputs;
    } Trellis;

    Trellis poly2trellis(const arma::imat &constraints,
        const arma::imat &codegen);
    Trellis poly2trellis(const arma::imat &constraints, 
        const arma::imat &codegen,
        const arma::imat &FeedbackConnection);
    
    arma::imat getHardDecisions(const arma::mat &x);

};

#endif
