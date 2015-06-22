// ballet includes
#include "viterbi.h"
#include "viterbi_p.h"

// system includes
#include <limits>

using namespace arma;

namespace ballet
{

    ViterbiDecoderPrivate::ViterbiDecoderPrivate()
        : balletObjectPrivate()
    {
    }

    ViterbiDecoderPrivate::~ViterbiDecoderPrivate()
    {
    }

    void ViterbiDecoderPrivate::lock()
    {
        BALLET_Q(ViterbiDecoder);

        // lock TrellisStructure property
        TrellisStructure = q->TrellisStructure;

        // lock InputFormat property
        InputFormat = q->InputFormat;

        // lock SoftInputWordLength
        SoftInputWordLength = q->SoftInputWordLength;

        // lock OutputFormat property
        OutputFormat = q->OutputFormat;

        // number of bits per symbol
        double kdbl = std::log((double)TrellisStructure.numInputSymbols) / std::log(2.0);
        k = static_cast<size_t>(kdbl);
        double ndbl = std::log((double)TrellisStructure.numOutputSymbols) / std::log(2.0);
        n = static_cast<size_t>(ndbl);

        // number of tail bits
        double numTails_ = log((double)TrellisStructure.numStates)/log(2.0);
        numTails = static_cast<int>(numTails_+0.5);

        // calculate bit output values
        binout.set_size(TrellisStructure.numOutputSymbols,n);
        if (InputFormat.compare("Unquantized") == 0)
        {
            for (size_t i=0; i<TrellisStructure.numOutputSymbols; i++)
            {
                int out_symbol = i;
                for (size_t j=0; j<n; j++)
                {
                    binout(i,j) = -2.0*static_cast<double>(out_symbol&1) + 1.0;
                    out_symbol = out_symbol >> 1;
                }
            }
        }
        else if (InputFormat.compare("Hard") == 0)
        {
            for (size_t i=0; i<TrellisStructure.numOutputSymbols; i++)
            {
                int out_symbol = i;
                for (size_t j=0; j<n; j++)
                {
                    binout(i,j) = static_cast<double>(out_symbol&1);
                    out_symbol = out_symbol >> 1;
                }
            }
        }
        else if (InputFormat.compare("Soft") == 0)
        {

            soft_zero = 0.0;
            soft_one = static_cast<double>((2<<SoftInputWordLength)-1);

            for (size_t i=0; i<TrellisStructure.numOutputSymbols; i++)
            {
                int out_symbol = i;
                for (size_t j=0; j<n; j++)
                {
                    binout(i,j) = soft_one * static_cast<double>(out_symbol&1);
                    out_symbol = out_symbol >> 1;
                }
            }

        }

        // set property_locked flag
        property_locked = true;
        
    }

    imat ViterbiDecoderPrivate::decode(const mat &x)
    {

        double inf = std::numeric_limits<double>::max();

        // number of trellis states
        size_t numStates = TrellisStructure.numStates;

        // coded frame length
        size_t cdeLen = x.n_elem;

        // number of trellis sections
        size_t L = cdeLen/n;

        // uncoded frame length
        size_t frmLen = (L-numTails) * k;

        // number number of trellis sections
        size_t NSIZ = numStates * (L+1);

        // [branch metric unit]
        mat branchmetric;
        if (InputFormat.compare("Unquantized") == 0)
        {
            branchmetric = branch_metric_unq(x);
        }
        else if (InputFormat.compare("Hard") == 0)
        {
            branchmetric = branch_metric_har(x);
        }
        else if (InputFormat.compare("Soft") == 0)
        {
            branchmetric = branch_metric_sof(x);
        }

        // path metric unit memory allocation
        mat phi(numStates,L+1);
        imat xi(numStates,L);
        imat lambda(numStates,L);

        // [path metric unit]
        {

            // [initialization]
            phi.fill(-inf);
            phi(0) = 0;
            lambda.zeros();

            for (size_t i=0; i<L; i++)
            {
                for (size_t s=0; s<numStates; s++)
                {
                    for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                    {

                        // calculate path metric
                        int nxt = TrellisStructure.nextStates(s,j);
                        int idx = TrellisStructure.outputs(s,j);
                        double mtr = branchmetric(idx,i) + phi(s,i);

                        // update and store path metric
                        if (mtr > phi(nxt,i+1))
                        {
                            phi(nxt,i+1) = mtr;
                            lambda(nxt,i) = j;
                            xi(nxt,i) = s;
                        }

                    }

                }

            }

        }

        // allocate traceback array
        imat symbols001(L,1);

        // [traceback]
        {

            // [termination state]
            int curstate = 0;

            for (size_t i=L; i>0; i--)
            {

                // most likely input symbol
                symbols001(i-1) = lambda(curstate,i-1);

                // traceback to previous trellis stage
                curstate = xi(curstate,i-1);

            }

        }

        // declare return value
        imat y;

        // format return value
        if (OutputFormat.compare("Binary") == 0)
        {

            // resize return matrix
            y.set_size(frmLen,1);
            
            for (size_t i=1; i<=L; ++i)
            {
    
                // decoded symbol value
                int input_symbol = symbols001(i-1);
   
                // convert symbol to binary vector 
                for (size_t kk=0; kk<k; ++kk)
                {
                    y(i*k-kk-1) = static_cast<int>(input_symbol&1);
                    input_symbol = input_symbol >> 1;
                }

            }

        } 
        else if (OutputFormat.compare("Integer") == 0)
        {

            // copy decoded symbols to output array
            y = symbols001(span(0,frmLen/k-1),0);

        }

        return y;

    }
    
    mat ViterbiDecoderPrivate::branch_metric_unq(const mat &x)
    {

        // number of trellis stages
        size_t len = x.n_elem / n;

        // allocate branch metric array
        mat y(TrellisStructure.numOutputSymbols,len);

        // [branch metric]
        for (size_t i=0; i<len; i++)
        {

            // get submatrix span
            mat x_ = flipud(x(span(i*n,(i+1)*n-1),0));

            // compute distance metric for all input combinations
            for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
            {

                // compute branch metric
                mat tmp = binout.row(j) * x_;

                // add branch metric to result
                y(j,i) = tmp(0);

            }

        }

        return y;

    }

    mat ViterbiDecoderPrivate::branch_metric_har(const mat &x)
    {

        // calculate number of trellis stages
        size_t len = x.n_elem / n;

        // allocate branch metric array
        mat y(TrellisStructure.numOutputSymbols,len);

        // convert binout to imat type 
        // [for hard branch metric computation]
        imat binout_ = conv_to<imat>::from(binout);

        // [branch metric]
        for (size_t i=0; i<len; i++)
        {

            // get submatrix span
            imat x_ = conv_to<imat>::from(
                flipud(x(span(i*n,(i+1)*n-1),0)));

            // compute distance metric for all input combinations
            for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
            {

                // compute exclusive or [XOR]
                imat bmet = binout_.row(j) + x_;

                // compute hard branch metric
                double tmp = 0.0;
                for (size_t nn=0; nn<n; ++nn)
                {
                    tmp += static_cast<double>(bmet(nn)%2);
                }

                y(j,i) = tmp;

            }

        }

        return y;

    }

    mat ViterbiDecoderPrivate::branch_metric_sof(const mat &x)
    {

        // calculate number of trellis stages
        size_t len = x.n_elem / n;

        // allocate branch metric array
        mat y(TrellisStructure.numOutputSymbols,len);

        // [branch metric]
        for (size_t i=0; i<len; i++)
        {

            // get buffer of bits for branch metric computation
            mat buf = trans(flipud(x(span(i*n,(i+1)*n-1),0)));

            // compute distance metric for all input combinations
            for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
            {

                // compute soft metric
                double bmet = n*soft_one - sum(buf-binout.row(j));
    
                // add branch metric to result
                y(j,i) = bmet;

            }

        }

        return y;

    }


    ViterbiDecoder::ViterbiDecoder()
        : balletObject(*new ViterbiDecoderPrivate)
    {

        // default trellis
        imat constraints = "7";
        imat codegen = "177 133;";
        TrellisStructure = poly2trellis(constraints,codegen);

        // default InputFormat
        InputFormat = "Unquantized";

        // default OutputFormat
        OutputFormat = "Integer";

    }

    imat ViterbiDecoder::decode(const mat &x)
    {

        BALLET_D(ViterbiDecoder);

        // lock decoder object
        if (!isLocked()) { d->lock(); }

        return d->decode(x);

    }

};
