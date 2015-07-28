#include "oct2dec.h"

using namespace arma;

namespace ballet
{

    imat oct2dec(const imat c)
    {

        // create decimal matrix d with same dimensions
        // as input octal matrix c
        imat d(c.n_rows,c.n_cols,fill::zeros);

        // iterate over input, output arrays
        {
            imat::const_iterator c_iter = c.begin();
            imat::iterator d_iter = d.begin();
            for ( ; c_iter != c.end(); ++c_iter, ++d_iter)
            {

                // get element value [octal]
                int oct_val = *c_iter;

                // initialize base
                int base = 1;

                // initialize result value [decimal]
                int dec_val = 0;
                
                while (oct_val != 0)
                {
            
                    // least significant digit
                    int oct_digit = oct_val%10;

                    // update decimal value
                    dec_val += oct_digit * base;

                    // get next octal digit
                    oct_val = oct_val/10;

                    // update base
                    base <<= 3;

                }

                // update result
                *d_iter = dec_val;

            }

        }

        return d;

    }

};
