#ifndef __ballet_lva_p_h
#define __ballet_lva_p_h

#include "lva.h"
#include "private/object_p.h"

namespace ballet
{

    namespace ListViterbiPrivate
    {

        typedef struct
        {
            int state;
            double phi;
            double chi;
            int lmbda;
        } path_t;

        typedef struct
        {
            splib::Trellis TrellisStructure;
            std::string OutputFormat;
            int L;
            size_t M;
            size_t k;
            int numTails;
            splib::vec branchmetric;
            std::vector<std::vector<path_t> > path_matrix;
        } ObjectData;

        typedef struct
        {
            splib::ivec merge_count;
            splib::imat merge_matrix;
            splib::vec phi;
            splib::vec chi;
            splib::ivec xi;
            splib::ivec lmbda;
        } SerialObjectData;

    };

    struct VDExtra
    {

        size_t k;           // number of bits per input symbol
        size_t n;           // number of bits per output symbol
        
        double soft_zero;   // most confident logical zero
        double soft_one;    // most confident logical one

        splib::vec binout;  // output codeword bits

    };

    class ListViterbiResultPrivate : public balletObjectPrivate
    {

        BALLET_DECLARE_PUBLIC(ListViterbiResult)

    public:
        ListViterbiResultPrivate();
        ~ListViterbiResultPrivate();

    public:
        splib::ivec next();
        splib::ivec operator()(const int n);

    public:
        int currentIndex();
        int count();

    public:
        ListViterbiPrivate::ObjectData * objectData;
        ListViterbiPrivate::SerialObjectData * serialObjectData;

    public:
        size_t num_paths;
        size_t current;

    };

    class ListViterbiDecoderPrivate : public balletObjectPrivate
    {

        BALLET_DECLARE_PUBLIC(ListViterbiDecoder)

    public:
        ListViterbiDecoderPrivate();
        ~ListViterbiDecoderPrivate();

    public:
        splib::Trellis TrellisStructure;
        std::string InputFormat;
        size_t SoftInputWordLength;
        std::string TerminationMethod;
        std::string OutputFormat;
        std::string Algorithm;
        int L;

    public:
        int numTails;

    public:
        ListViterbiResult decode(const splib::fvec &x);
        void lock();

    public:
        void createExtra();
        void deleteExtra();

    public:
        VDExtra *extra;

    protected:
        ListViterbiResult decode_parallel(const splib::fvec &x);
        ListViterbiResult decode_serial(const splib::fvec &x);

    protected:
        splib::vec branch_metric_unq(const splib::fvec &x);
        splib::vec branch_metric_har(const splib::fvec &x);
        splib::vec branch_metric_sof(const splib::fvec &x);

    };

};

#endif
