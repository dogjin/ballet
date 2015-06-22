#include "logical.h"

using namespace arma;

namespace ballet
{

    int all(const umat & A)
    {

        // initialize function result
        int result = 1;

        // iterate over array
        for (umat::const_iterator i=A.begin(); i!=A.end(); i++)
        {
            if (*i == 0)
            {
                result = 0;
                break;
            }
        }

        return result;

    }

};
