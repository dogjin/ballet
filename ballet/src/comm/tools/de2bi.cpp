#include "de2bi.h"
#include "de2bi_p.h"

using namespace arma;

namespace ballet
{

    imat de2biPrivate::de2bi(const imat &d, const int n)
    {

        // create return matrix
        imat b(d.n_elem,n,fill::zeros);

        // iterate over input, output arrays
        int current_row = 0;
        for (imat::iterator d_iter = d.begin();
             d_iter!=d.end(); ++d_iter)
        {

            // get element value [decimal]
            int dec_val = *d_iter;

            // row iterator for i:th binary value
            mat::row_iterator b_iter = b.begin_row(current_row);

            while (dec_val != 0)
            {

                // least significant digit
                *b_iter = dec_val%2;

                // get next binary digit
                dec_val >>= 1;

                // update row iterator, next column
                ++b_iter;

            }

            // update current row
            ++current_row;

        }

        return b;

    }

    imat de2bi(const int _d)
    {

        imat d = _d;
        return de2bi(d);

    }

    imat de2bi(const int _d, const int n)
    {
        imat d = _d;
        retrun de2bi(d,n);
    }

    imat de2bi(const imat d)
    {

        // get largest decimal value
        int max_val = max(max(d));

        // calculate number of bits needed to convert
        // decimal matrix d to binary matrix b
        int ntmp;
        if (max_val != 0)
        {

            double mv_dbl = static_cast<double>(max_val);
            double e = std::floor(std::log2(mv_dbl)) + 1.0;
            ntmp = static_cast<int>(e);
    
        }
        else
        {
            ntmp = 1;
        }

        return de2biPrivate::de2bi(d,n);

    }

    imat de2bi(const imat d, const int n)
    {

        return de2biPrivate(d,n);


    }

};
