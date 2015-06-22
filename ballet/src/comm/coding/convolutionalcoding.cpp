#include "convolutionalcoding.h"
#include <ballet/comm/de2bi.h>
#include <ballet/comm/bi2de.h>
#include <ballet/comm/oct2dec.h>

using namespace arma;

namespace ballet
{

    Trellis poly2trellis(const imat &constraints, const imat &codegen)
    {
        return Trellis();
    }

    Trellis poly2trellis(const imat &constraints, const imat &codegen, const imat &FeedbackConnection)
    {

        // get k, n values for rate k/n convolutional code
        int k = constraints.n_elem;
        int n = codegen.n_cols;

        // constraint length
        int K = sum(vectorise(constraints));

        // number of memory elements
        int L = K - k;

        // number of states, symbols and transistions
        int NumberStates = 1 << L;
        int NumberSymbols = k << 1;
        int number_transitions = 1 << K;

        /*
         * Build the encoder trellis with the following information:
         *
         * - Initial state number, init (element of {0,2,...,M^L-1})
         * - Final state number, final (element of {0,2,...,M^L-1})
         * - new symbol, ns
         * - codeword value, cword
         *
         * Enumerate all possible transitions and assign corresponding input symbols
         */
        imat init(number_transitions,1);
        imat final(number_transitions,1);
        imat ns(number_transitions,1);
        imat cword(number_transitions,1);

        for (int ii=0; ii<number_transitions; ii++)
        {

            // convert transition to input / register
            // value binary representation
            imat btemp = de2bi(ii,K);

            // initialize internals
            int indx = 0;
            imat bvec_temp(n,1,fill::zeros);
            imat init_tmp;
            imat final_tmp;
            imat ns_tmp;

            /*
             * Generate codeword output for each set of input registers.
             * There are two different ways to generate the state vectors. Lin and Costello start
             * with the 1st set of registers as the most significant bits and decrease in a 
             * downward fashion. However, we can also start with the last set of registers 
             * as the most significant bits and decrease upwards. Both methods yield the same results
             * in the final encoding and decoding.
             */
            for (int kk=0; kk<k; kk++)
            {

                // resiter (kk)
                int Ktemp = constraints(kk,0);
                int Ltemp = Ktemp - 1;
                imat bits = btemp(span(indx,indx+Ktemp-1),0);
                init_tmp = join_cols(init_tmp,bits(span(1,Ktemp-1),0));
                ns_tmp = join_cols(ns_tmp,bits(span(0,0),0));

                // apply FeedbackConnection(kk)
	            imat fb_temp = FeedbackConnection;	
                imat fb_dval = oct2dec(fb_temp);
                fb_temp = de2bi(fb_dval,Ktemp);
                imat bitval = fb_temp * bits;
                bits = join_cols(bitval,btemp(span(indx+1,indx+Ktemp-1),0));

                // update trellis transition values with
                // results from shift register (kk)
                final_tmp = join_cols(final_tmp,bits(span(0,Ltemp-1),0));
                imat tempv = oct2dec(codegen.row(kk));
                bvec_temp = bvec_temp + de2bi(tempv)*bits;
                indx = Ktemp;

            }

            // store results of state transition
            init(span(ii,ii),0) = bi2de(init_tmp);
            final(span(ii,ii),0) = bi2de(final_tmp);
            ns(span(ii,ii),0) = bi2de(ns_tmp);
            cword(span(ii,ii),0) = bi2de(bvec_temp);

        }

        // map transition representation to
        // trellis state representation
        imat StateMap(NumberStates,k<<1);
        imat EncodeMap(NumberStates,k<<1);
        for (int ii=0; ii<number_transitions; ii++)
        {
            StateMap(init(ii),ns(ii)) = final(ii);
            EncodeMap(init(ii),ns(ii)) = cword(ii);
        }

        // convert to trellis representation
        Trellis trellis;
        trellis.nextStates = StateMap;
        trellis.outputs = EncodeMap;
        trellis.numInputSymbols = NumberSymbols;
        trellis.numOutputSymbols = static_cast<int>(std::pow((float)2,(float)n));
        trellis.numStates = NumberStates;

        return trellis;

    }

    imat getHardDecisions(const mat &x)
    {

        // unquantized vector length
        int L = x.n_elem;

        // declare quanized vector
        imat y(L,1);

        // convert unquantized vector
        // to quantized bits
        for (size_t i=0; i<L; i++)
        {
            if (x(i) > 0)
                y(i) = 1;
            else
                y(i) = 0;
        }

        return y;

    }

};
