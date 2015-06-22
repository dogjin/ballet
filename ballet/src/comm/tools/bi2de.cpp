#include "bi2de.h"

using namespace arma;

namespace ballet
{

    imat bi2de(const imat b)
    {

        // create return vector
        imat d(b.n_rows,1);

        int row_number = 0;
        for (size_t i=0; i<b.n_rows; ++i)
        {

            // initialize base
            int base = 1;

            // initialize result value [decimal]
            int dec_val = 0;

            // row iterator for i:th binary value
            imat::row_iterator b_iter = b.begin_row(i);
           
            // convert binary value to decimal number 
            while (b_iter != b.end_row(i))
            {

                if (*b_iter == 0) { base <<= 1; continue; }
                dec_val += base;
                
            }

            // update result
            d(i) = dec_val;

        }

        return d;

    }

};
