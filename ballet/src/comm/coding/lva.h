#ifndef __ballet_comm_lva_h
#define __ballet_comm_lva_h

#include <vector>
#include <splib/vec.h>
#include <splib/commfunc.h>
#include <ballet/object.h>

typedef signed char qint8;

namespace ballet
{

    class ListViterbiDecoder;
    class ListViterbiDecoderPrivate;
    class ListViterbiResultPrivate;

    class ListViterbiResult : public balletObject
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(ListViterbiResult)
        friend class ListViterbiDecoder;
        friend class ListViterbiDecoderPrivate;
        /// @endcond

    public:
        ListViterbiResult();
        ListViterbiResult(const ListViterbiResult &other);

    public:
        ListViterbiResult& operator=(const ListViterbiResult &other);

    protected:
        ListViterbiResult(size_t n);

    public:
        
        splib::ivec next();
        splib::ivec operator()(const int n);

    protected:
        std::vector<splib::ivec> p;

    };

    class ListViterbiDecoder : public balletObject
    {

        /// @cond internal
        BALLET_DECLARE_PRIVATE(ListViterbiDecoder)
        /// @endcond

    public:
        ListViterbiDecoder();

    public:
        splib::Trellis TrellisStructure;
        std::string InputFormat;
        size_t SoftInputWordLength;
        std::string TerminationMethod;
        std::string OutputFormat;
        std::string Algorithm;
        int L;

    public:
        ListViterbiResult decode(const splib::fvec &x);

    };

};

#endif
