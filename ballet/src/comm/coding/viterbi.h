#ifndef __ballet_comm_viterbi_h
#define __ballet_comm_viterbi_h

// ballet includes
#include "ballet/object.h"
#include "convolutionalcoding.h"

// external includes
#include <armadillo>

namespace ballet
{

    class ViterbiDecoderPrivate;

    class ViterbiDecoder : public balletObject
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(ViterbiDecoder)
        /// @endcond

    public:
        ViterbiDecoder();

    public:
        Trellis TrellisStructure;
        std::string InputFormat;
        size_t SoftInputWordLength;
        std::string OutputFormat;

    public:
        arma::imat decode(const arma::mat &x);

    };

};

#endif
