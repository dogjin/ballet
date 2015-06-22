#include <armadillo>
#include "random.h"

using namespace arma;

namespace ballet
{

    ivec randperm(const int n)
    {

        // generate uniform random number sequence
        vec x = randu<vec>(n);

        // sort uniform random number sequence
        // to get random permutation vector
        uvec y = sort_index(x);

        return conv_to<ivec>::from(y);

    }

};
