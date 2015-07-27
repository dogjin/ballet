#define TWOPI 6.283185307
#define Billion 1.e9

// ballet includes
#include "appdecoder.h"
#include "private/object_p.h"

// system includes
#include <vector>
#include <limits>

#define Poly(x) (0.6949317915+x*(-0.5196305002+x*(0.1617251897+x*(-0.0257407636+x*(0.0020686227+x*(-0.0000665807))))))
#define correction(x) (x<8?Poly((x)):3.35406e-4)

using namespace arma;

namespace ballet
{

    typedef struct
    {
        size_t count;
        std::vector<int> inpt;
        std::vector<int> outp;
        std::vector<int> prev;
    } pstate_t;

    typedef struct
    {
        size_t k;               // number of bits per input symbol
        size_t n;               // number of bits per output symbol
        mat bininp;             // vector of input bit codewords
        struct next_t
        {
            mat binout;         // vector of output bit codewords
        } next;
    } enc_t;

    class APPDecoderWorkerBase
    {
    public:
        APPDecoderWorkerBase()  {}
        ~APPDecoderWorkerBase() {}

    public:
        Trellis TrellisStructure;
        std::string Algorithm;
        std::string TerminationMethod;
        size_t NumScalingBits;
        double ScalingFactor;

    public:
        enc_t enc;
        std::vector<pstate_t> pStateMap;

    public:
        virtual void decode(const mat &lu, const mat &lc, mat *lcd, mat *lud) {}

    };

    class APPDecoderWorker_MaxStar
        : public APPDecoderWorkerBase
    {

    public:
        void decode(const mat &lu, const mat &lc, mat *lcd, mat *lud);

    };

    class APPDecoderWorker_LogMap 
        : public APPDecoderWorkerBase
    {

    public:
        void decode(const mat &lu, const mat &lc, mat *lcd, mat *lud);

    };

    class APPDecoderWorker_Max 
        : public APPDecoderWorkerBase
    {

    public:
        void decode(const mat &lu, const mat &lc, mat *lcd, mat *lud);

    };

    /*
      Convenience class for APP Decoder
    */
    class APPDecoderPrivate : public balletObjectPrivate
    {

        BALLET_DECLARE_PUBLIC(APPDecoder)

    public:
        APPDecoderPrivate();
        ~APPDecoderPrivate();

    public:
        std::string Algorithm;

    public:
        void decode(const mat &lu, const mat &lc, mat *lcd, mat *lud);

    public:
        void lock();

    public:

        /*!
          Computes Max* Log-MAP simplification
        */
        static float maxStar(const float &x, const float &y)
        {

            return std::max(x,y) + correction(fabs(x-y));

        }

        /*!
          Computes Log-MAP
        */
        static float logMAP(const float &x, const float &y)
        {
            float out = std::max(x,y) + std::log(1.0+std::exp(-std::abs(x-y)));
            return out;
        }

        /*!
          Computes Max-Log-MAP
        */
        static float max_(const float &x, const float &y)
        {
            return std::max(x,y);
        }

    private:
        static float maxstar_table[176];

    private:

       APPDecoderWorkerBase * worker;

    };

    APPDecoderPrivate::APPDecoderPrivate()
        : balletObjectPrivate()
          , worker(0)
    {}

    APPDecoderPrivate::~APPDecoderPrivate()
    {
        if (worker)
        {
            delete worker;
            worker = 0;
        }
    }

    float APPDecoderPrivate::maxstar_table[] = {
        0.69314718055994529, 0.673347167228034, 0.65394696731758994, 0.63494610159561349, 
        0.61634377304073962, 0.59813886938159178, 0.58032996662642578, 0.5629153335603464, 
        0.54589293718007537, 0.52926044903028424, 0.5130152523999526, 0.49715445033210992, 
        0.48167487439574341, 0.46657309416461801, 0.45184542734430633, 0.43748795048588568, 
        0.42349651022253426, 0.4098667349636621, 0.39659404698022444, 0.38367367481449383, 
        0.37110066594777769, 0.35886989966032323, 0.3469761000189524, 0.33541384892973053, 
        0.32417759919518879, 0.31326168751822286, 0.30266034739773889, 0.29236772186435833, 
        0.28237787600797598, 0.27268480925263938, 0.26328246733803129, 0.25416475397074761, 
        0.24532554211251709, 0.23675868487646651, 0.22845802600646797, 0.22041740991845082, 
        0.21263069128632331, 0.20509174415876136, 0.19779447059659636, 0.19073280882382182, 
        0.18390074088833888, 0.17729229983146019, 0.17090157636787057, 0.16472272508020849, 
        0.15874997013467182, 0.15297761052607406, 0.14740002486257028, 0.14201167570185888, 
        0.13680711345203822, 0.13178097985146936, 0.1269280110429726, 0.12224304025848917, 
        0.11772100013096001, 0.11335692465064118, 0.10914595078339807, 0.10508331976869593, 
        0.10116437811507242, 0.097384578310816469, 0.093739479267430328, 0.0902247465132089, 
        0.086836152153949631, 0.083569574617418818, 0.080420998197756693, 0.077386512415507938, 
        0.074462311208430443, 0.071644691967669719, 0.068930054433295349, 0.066314899462582039, 
        0.063795827683805567, 0.061369538047684018, 0.059032826287971366, 0.056782583302082891, 
        0.054615793462002286, 0.052529532865117079, 0.050520967534021702, 0.048587351573741958, 
        0.046726025294271521, 0.04493441330574701, 0.043210022593073764, 0.041550440576283008, 
        0.039953333162430306, 0.038416442794361086, 0.036937586501232821, 0.035514653955253266, 
        0.034145605538694994, 0.032828470424865447, 0.031561344676348559, 0.030342389363505955, 
        0.029169828705895982, 0.028041948238979979, 0.02695709300820814, 0.025913665792307073, 
        0.024910125357366336, 0.023944984743078702, 0.023016809582299336, 0.022124216454879247, 
        0.021265871276566935, 0.020440487723596214, 0.019646825693436724, 0.018883689802042317, 
        0.018149927917809779, 0.017444429732341199, 0.016766125368008675, 0.016113984022215113, 
        0.015487012648170265, 0.014884254671918118, 0.014304788745287698, 0.013747727534377217, 
        0.013212216543127738, 0.012697432971496305, 0.012202584607696109, 0.011726908753935301, 
        0.011269671185057697, 0.010830165139457164, 0.010407710341623785, 0.01000165205565178, 
        0.0096113601690349156, 0.0092362283060555637, 0.0088756729700722216, 0.0085291327139978453, 
        0.0081960673382677641, 0.0078759571155826037, 0.0075683020417261944, 0.0072726211117517415, 
        0.0069884516208369744, 0.0067153484891179669, 0.006452883609813876, 0.0062006452199646431, 
        0.0059582372931190029, 0.0057252789533070959, 0.0055014039096574737, 0.0052862599110215227, 
        0.0050795082199807731, 0.0048808231056280821, 0.0046898913545248347, 0.0045064117992494046, 
        0.0043300948639671951, 0.0041606621264625599, 0.0039978458960905749, 0.0038413888071198794, 
        0.0036910434269464489, 0.0035465718786807029, 0.0034077454776149101, 0.0032743443810994677, 
        0.0031461572513634024, 0.0030229809308315535, 0.0029046201295047304, 0.0027908871239777167, 
        0.0026816014676889263, 0.0025765897120008882, 0.0024756851377303571, 0.0023787274967537238, 
        0.0022855627633260227, 0.0021960428947675734, 0.0021100256011755743, 0.0020273741238382385, 
        0.0019479570220328406, 0.0018716479679019994, 0.0017983255491144409, 0.0017278730790231504, 
        0.0016601784140455901, 0.0015951337780006952, 0.0015326355931442586, 0.0014725843176540952, 
        0.001414884289328097, 0.0013594435752599973, 0.0013061738272732621, 0.0012549901428946441, 
        0.001205810931664377, 0.001158557786576929, 0.0011131553604645742, 0.0010695312471351985, 
        0.0010276158670836635, 0.00098734235760847559, 0.00094864646716139598, 0.00091146645377420147
    };

    /*!
      Lock nontunable properties and input specifications
    */
    void APPDecoderPrivate::lock()
    {

        if (isLocked()) { return; }

        BALLET_Q(APPDecoder);

        // destroy worker instance [if applicable]
        if (worker) { delete worker; worker = 0; }

        // lock algorithm
        Algorithm = q->Algorithm;
        if (Algorithm.compare("True APP") == 0)
        {
            worker = new APPDecoderWorker_LogMap;
            worker->Algorithm = Algorithm;
        }
        else if (Algorithm.compare("Max*") == 0)
        {
            worker = new APPDecoderWorker_MaxStar;
            worker->Algorithm = Algorithm;
        }
        else if (Algorithm.compare("Max") == 0)
        {
            worker = new APPDecoderWorker_Max;
            worker->Algorithm = Algorithm;
        }
        else
        {
            printf("Unknown Decoding Algorithm property in APPDecoder, defaulting to <Max*>\n");
            Algorithm = "Max*";
            worker = new APPDecoderWorker_MaxStar;
            worker->Algorithm = Algorithm;
        }

        // lock trellis
        worker->TrellisStructure = q->TrellisStructure;

        // lock encoder
        {

            // number of bits per input symbol
            double kdbl = std::log((double)worker->TrellisStructure.numInputSymbols) / std::log(2.0);
            worker->enc.k = static_cast<size_t>(kdbl);

            // number of bits per output symbol
            double ndbl = std::log((double)worker->TrellisStructure.numOutputSymbols) / std::log(2.0);
            worker->enc.n = static_cast<size_t>(ndbl);

            // create vector of bit input values
            worker->enc.bininp.set_size(worker->TrellisStructure.numInputSymbols,worker->enc.k);

            // calculate bit input values
            for (size_t i=0; i<worker->TrellisStructure.numInputSymbols; i++)
            {
                int input_symbol = static_cast<int>(i);
                for (size_t j=0; j<worker->enc.k; j++)
                {
                    worker->enc.bininp(i,j) = 2.0*static_cast<double>(input_symbol&1) - 1.0;
                    input_symbol = input_symbol >> 1;
                }
            }

            // create vector of bit output values
            worker->enc.next.binout.set_size(worker->TrellisStructure.numOutputSymbols,worker->enc.n);

            // calculate bit output values
            for (size_t i=0; i<worker->TrellisStructure.numOutputSymbols; i++)
            {
                int out_symbol = static_cast<int>(i);
                for (size_t j=0; j<worker->enc.n; j++)
                {
                    worker->enc.next.binout(i,j) = 2.0*static_cast<double>(out_symbol&1) - 1.0;
                    out_symbol = out_symbol >> 1;
                }
            }

        }

        // lock TerminationMethod
        worker->TerminationMethod = q->TerminationMethod;
        if ( worker->TerminationMethod.compare("Terminated") && worker->TerminationMethod.compare("Truncated"))
        {
            printf("TerminationMethod property in APPDecoder not one of Terminated|Truncated, defaulting to Terminated\n");
            worker->TerminationMethod = "Terminated";
        }

        // lock NumScalingBits
        worker->NumScalingBits = q->NumScalingBits;
        if (q->NumScalingBits < 0 || q->NumScalingBits > 8)
        {
            printf("NumScalingBits property in APPDecoder outside range, defaulting to 3\n");
            worker->NumScalingBits = 3;
        }
        worker->ScalingFactor = static_cast<float>(1<<worker->NumScalingBits);

        // initialize previous state map
        worker->pStateMap.resize(worker->TrellisStructure.numStates);
        for (size_t s=0; s<worker->TrellisStructure.numStates; s++)
        {
            worker->pStateMap[s].count = 0;
            worker->pStateMap[s].inpt.clear();
            worker->pStateMap[s].outp.clear();
            worker->pStateMap[s].prev.clear();
        }

        // compute previous state map
        {

            for (size_t s=0; s<worker->TrellisStructure.numStates; s++)
            {
                for (size_t j=0; j<worker->TrellisStructure.numInputSymbols; j++)
                {

                    // calculate path metric
                    int nxt = worker->TrellisStructure.nextStates(s,j);
                    int idx = worker->TrellisStructure.outputs(s,j);

                    pstate_t & pState = worker->pStateMap[nxt];
                    pState.prev.push_back(s);
                    pState.inpt.push_back(j);
                    pState.outp.push_back(idx);
                    pState.count++;

                }
            }

        }

        // set lock flag
        property_locked = true;

    }

    void APPDecoderPrivate::decode(const mat &Lu, const mat &Lc, mat *lcd, mat *lud)
    {

        return worker->decode(Lu,Lc,lcd,lud);

    }

    // ***************************** Algorithm = LogMap ************************************************ //
    void APPDecoderWorker_LogMap::decode(const mat &Lu, const mat &Lc, mat *lcd, mat *lud)
    {

        // number of trellis states
        size_t numStates = TrellisStructure.numStates;

        // maximum floating point value
        double inf = std::numeric_limits<double>::max();

        // calculate output, traceback and traceforward matrix sizes
        size_t len = Lc.n_elem / enc.n;

        // traceforward path metric matrix
        mat ak (numStates,len+1);
        ak.fill(-inf); ak(0) = 0;

        // traceback path metric matrix
        mat bk (numStates,len+1);
        bk.fill(-inf);
        if (TerminationMethod.compare("Terminated") == 0)
        {
            bk(0,len) = 0;
        }
        else if (TerminationMethod.compare("Truncated") == 0)
        {
            bk(span::all,len) = zeros<mat>(numStates,1);
        }

        // time-dependent scaling factor to prevent
        // excessive growth of the numerical values of
        // alpha and beta matricies
        vec denom(len); denom.fill(-inf);

        // branch transition likelihood matricies
        mat gammau(TrellisStructure.numInputSymbols,len);
        mat gammac(TrellisStructure.numOutputSymbols,len);

        // calculate uncoded transition probability values
        {

            for (size_t i=0; i<len; i++)
            {
            
                // compute distance metric for all input combinations
                for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                {

                    mat tmp = enc.bininp.row(j) * flipud(Lu(span(i*enc.k,(i+1)*enc.k),0));
                    gammau(j,i) = 0.5 * tmp(0);

                }

            }

        }

        // calculate coded transition probability values
        {

            for (size_t i=0; i<len; i++)
            {
            
                // compute distance metric for all input combinations
                for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
                {

                    mat tmp = enc.next.binout.row(j) *
                        flipud(Lc(span(i*enc.n,(i+1)*enc.n-1),
                        0));
                    gammac(j,i) = 0.5 * tmp(0);

                }

            }

        }

        // calculate alpha values [forward recursion]
        {

            for (size_t i=0; i<len; i++)
            {

                for (size_t s=0; s<numStates; s++)
                {

                    pstate_t & pState = pStateMap[s];

                    if (pState.count == 0)
                    {
                        ak(numStates,i) = -inf;
                        continue;
                    }

                    ak(numStates,i) = gammau(pState.inpt[0],i) + gammac(pState.outp[0],i) + ak(pState.prev[0],i);

                    for (size_t j=1; j<pState.count; j++)
                    {
                        double mtr = gammau(pState.inpt[j],i) + gammac(pState.outp[j],i) + ak(pState.prev[j],i);
                        ak(s,i+1) = APPDecoderPrivate::logMAP(mtr,ak(s,i+1));
                    }

                    denom(i) = APPDecoderPrivate::logMAP(ak(s,i+1),denom(i));

                }

                // normalize alpha
                for (size_t s=0; s<numStates; s++)
                {
                    ak(s,i+1) -= denom(i);
                }

            }

        }

        // calculate beta values [backwards recursion]
        {

            for (size_t i=len; i>0; i--)
            {

                for (size_t s=0; s<numStates; s++)
                {

                    bk(s,i-1) = gammau(0,len-1) + gammac(TrellisStructure.outputs(s,0)) + bk(TrellisStructure.nextStates(s,0),len);

                    for (size_t j=1; j<TrellisStructure.numInputSymbols; j++)
                    {
                        double mtr = gammau(j,len-1) + gammac(TrellisStructure.outputs(s,j)) + bk(TrellisStructure.nextStates(s,j),len);
                        bk(s,i-1) = APPDecoderPrivate::logMAP(mtr,bk(s,i-1));
                    }

                }

                // normalize beta
                for (size_t s = 0; s<numStates; s++)
                {
                    bk(s,i-1) -= denom(i-1);
                }

            }

        }

        // pre-allocate unencoded bit probability
        // output vector L(u)
        lud->set_size(Lu.n_rows,Lu.n_cols);

        // update unencoded bit probability
        // log-likelihood ratio
        {

            mat & LUD = *lud;

            for (size_t i=0; i<len; i++)
            {

                for (size_t k=0; k<enc.k; k++)
                {
                    
                    float one_val = -inf;
                    float zero_val = -inf;

                    int mask = (1<<(enc.k-k-1));

                    for (size_t s=0; s<numStates; s++)
                    {

                        for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                        {

                            // calculate path metric
                            int nxt = TrellisStructure.nextStates(s,j);
                            int idx = TrellisStructure.outputs(s,j);
                            double metr = ak(s,i) + bk(nxt,i+1) + gammac(idx,i);

                            if (j&mask)     // binary one input
                                one_val = APPDecoderPrivate::logMAP(metr,one_val);
                            else
                                zero_val = APPDecoderPrivate::logMAP(metr,zero_val);

                        }

                    }

                    // update output metric
                    LUD(i*enc.k+k) = one_val - zero_val;

                }

            }

        }

        // update codeded bit probability
        // log-likelihood ratio
        if ( lcd != 0 )
        {

            // resize encoeded bit probability
            // return vector
            lcd->set_size(Lc.n_rows,Lc.n_cols);

            // convert object pointer to object reference
	        mat & LCD = *lcd;	

            for (size_t i=0; i<len; i++)
            {
                
                for (size_t n=0; n<enc.n; n++)
                {
                    
                    double one_val = -inf;
                    double zero_val = -inf;

                    for (size_t s=0; s<numStates; s++)
                    {

                        for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                        {

                            // calculate path metric
                            int nxt = TrellisStructure.nextStates(s,j);
                            int idx = TrellisStructure.outputs(s,j);
                            double metr = ak(s,i) + bk(nxt,i+1) + gammau(j,i);

                            if (idx&(1<<(enc.n-n-1)))       // binary one input
                                one_val = APPDecoderPrivate::logMAP(metr,one_val);
                            else                            // binary zero input
                                zero_val = APPDecoderPrivate::logMAP(metr,zero_val);

                        }

                    }

                    // update output metric
                    LCD(i*enc.n+n) = one_val - zero_val;

                }

            }

        }

    }

    // ***************************** Algorithm = Max* ************************************************ //
    void APPDecoderWorker_MaxStar::decode(const mat &Lu, const mat &Lc, mat *lcd, mat *lud)
    {

        // number of trellis states
        size_t numStates = TrellisStructure.numStates;

        // scale input to avoid losing precision
        mat * lu_scaled = 0;
        mat * lc_scaled = 0;
        {

            // create new matrix objects
            // for scaled values
            mat * lu_scaled = new mat;
            mat * lc_scaled = new mat;

            // convert object pointer to 
            // object reference
            mat & tmp001 = *lu_scaled;
            mat & tmp002 = *lc_scaled;

            // scale uncoded log-likelihood
            // probability bit vector
            tmp001 = Lu * ScalingFactor;

            // scale codeded log-likelihood
            // probability bit vector
            tmp002 = Lc * ScalingFactor;

        }

        // maximum floating point value
        double inf = std::numeric_limits<double>::max();

        // calculate output, traceback and traceforward matrix sizes
        size_t len = Lc.n_elem / enc.n;

        // traceforward path metric matrix
        mat ak (numStates,len+1);
        ak.fill(-inf); ak(0) = 0;

        // traceback path metric matrix
        mat bk (numStates,len+1);
        bk.fill(-inf);
        if (TerminationMethod.compare("Terminated") == 0)
        {
            bk(0,len) = 0;
        }
        else if (TerminationMethod.compare("Truncated") == 0)
        {
            bk(span::all,len) = zeros<mat>(numStates,1);
        }

        // time-dependent scaling factor to prevent
        // excessive growth of the numerical values of
        // alpha and beta matricies
        vec denom(len); denom.fill(-inf);

        // branch transition likelihood matricies
        mat gammau(TrellisStructure.numInputSymbols,len);
        mat gammac(TrellisStructure.numOutputSymbols,len);

        // calculate uncoded transition probability values
        {

            // convert object pointer to
            // reference object
            mat & tmp001 = *lu_scaled;

            for (size_t i=0; i<len; i++)
            {
            
                // compute distance metric for all input combinations
                for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                {

                    mat tmp = enc.bininp.row(j) * flipud(tmp001(span(i*enc.k,(i+1)*enc.k),0));
                    gammau(j,i) = 0.5 * tmp(0);

                }

            }

        }

        // calculate coded transition probability values
        {

            // convert object pointer to
            // reference object
            mat & tmp001 = *lc_scaled;

            for (size_t i=0; i<len; i++)
            {
            
                // compute distance metric for all input combinations
                for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
                {

                    mat tmp = enc.next.binout.row(j) *
                        flipud(tmp001(span(i*enc.n,(i+1)*enc.n-1),
                        0));
                    gammac(j,i) = 0.5 * tmp(0);

                }

            }

        }

        // calculate alpha values [forward recursion]
        {

            for (size_t i=0; i<len; i++)
            {

                for (size_t s=0; s<numStates; s++)
                {

                    pstate_t & pState = pStateMap[s];

                    if (pState.count == 0)
                    {
                        ak(numStates,i) = -inf;
                        continue;
                    }

                    ak(numStates,i) = gammau(pState.inpt[0],i) + gammac(pState.outp[0],i) + ak(pState.prev[0],i);

                    for (size_t j=1; j<pState.count; j++)
                    {
                        double mtr = gammau(pState.inpt[j],i) + gammac(pState.outp[j],i) + ak(pState.prev[j],i);
                        ak(s,i+1) = APPDecoderPrivate::maxStar(mtr,ak(s,i+1));
                    }

                    denom(i) = APPDecoderPrivate::maxStar(ak(s,i+1),denom(i));

                }

                // normalize alpha
                for (size_t s=0; s<numStates; s++)
                {
                    ak(s,i+1) -= denom(i);
                }

            }

        }

        // calculate beta values [backwards recursion]
        {

            for (size_t i=len; i>0; i--)
            {

                for (size_t s=0; s<numStates; s++)
                {

                    bk(s,i-1) = gammau(0,len-1) + gammac(TrellisStructure.outputs(s,0)) + bk(TrellisStructure.nextStates(s,0),len);

                    for (size_t j=1; j<TrellisStructure.numInputSymbols; j++)
                    {
                        double mtr = gammau(j,len-1) + gammac(TrellisStructure.outputs(s,j)) + bk(TrellisStructure.nextStates(s,j),len);
                        bk(s,i-1) = APPDecoderPrivate::maxStar(mtr,bk(s,i-1));
                    }

                }

                // normalize beta
                for (size_t s = 0; s<numStates; s++)
                {
                    bk(s,i-1) -= denom(i-1);
                }

            }

        }

        // pre-allocate unencoded bit probability
        // output vector L(u)
        lud->set_size(Lu.n_rows,Lu.n_cols);

        // update unencoded bit probability
        // log-likelihood ratio
        {

            mat & LUD = *lud;

            for (size_t i=0; i<len; i++)
            {

                for (size_t k=0; k<enc.k; k++)
                {
                    
                    float one_val = -inf;
                    float zero_val = -inf;

                    int mask = (1<<(enc.k-k-1));

                    for (size_t s=0; s<numStates; s++)
                    {

                        for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                        {

                            // calculate path metric
                            int nxt = TrellisStructure.nextStates(s,j);
                            int idx = TrellisStructure.outputs(s,j);
                            double metr = ak(s,i) + bk(nxt,i+1) + gammac(idx,i);

                            if (j&mask)     // binary one input
                                one_val = APPDecoderPrivate::maxStar(metr,one_val);
                            else
                                zero_val = APPDecoderPrivate::maxStar(metr,zero_val);

                        }

                    }

                    // update output metric
                    LUD(i*enc.k+k) = one_val - zero_val;

                }

            }

        }

        // update codeded bit probability
        // log-likelihood ratio
        if ( lcd != 0 )
        {

            // resize encoeded bit probability
            // return vector
            lcd->set_size(Lc.n_rows,Lc.n_cols);

            // convert object pointer to object reference
	        mat & LCD = *lcd;	

            for (size_t i=0; i<len; i++)
            {
                
                for (size_t n=0; n<enc.n; n++)
                {
                    
                    double one_val = -inf;
                    double zero_val = -inf;

                    for (size_t s=0; s<numStates; s++)
                    {

                        for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                        {

                            // calculate path metric
                            int nxt = TrellisStructure.nextStates(s,j);
                            int idx = TrellisStructure.outputs(s,j);
                            double metr = ak(s,i) + bk(nxt,i+1) + gammau(j,i);

                            if (idx&(1<<(enc.n-n-1)))       // binary one input
                                one_val = APPDecoderPrivate::maxStar(metr,one_val);
                            else                            // binary zero input
                                zero_val = APPDecoderPrivate::maxStar(metr,zero_val);

                        }

                    }

                    // update output metric
                    LCD(i*enc.n+n) = one_val - zero_val;

                }

            }

        }


        {

            // scale uncoded log-likelihood bit
            // vector output by inverse of scaling factor
            *lud = *lud / ScalingFactor;

            // scale coded log-likelihood bit
            // vector output by inverse of scaling factor
            if ( lcd != 0 )
            {
                *lcd = *lcd / ScalingFactor;
            }

            delete lu_scaled;
            delete lc_scaled;

        }

    }


    // ***************************** Algorithm = Max ************************************************ //
    void APPDecoderWorker_Max::decode(const mat &Lu, const mat &Lc, mat *lcd, mat *lud)
    {

        // number of trellis states
        size_t numStates = TrellisStructure.numStates;

        // maximum floating point value
        double inf = std::numeric_limits<double>::max();

        // calculate output, traceback and traceforward matrix sizes
        size_t len = Lc.n_elem / enc.n;

        // traceforward path metric matrix
        mat ak (numStates,len+1);
        ak.fill(-inf); ak(0) = 0;

        // traceback path metric matrix
        mat bk (numStates,len+1);
        bk.fill(-inf);
        if (TerminationMethod.compare("Terminated") == 0)
        {
            bk(0,len) = 0;
        }
        else if (TerminationMethod.compare("Truncated") == 0)
        {
            bk(span::all,len) = zeros<mat>(numStates,1);
        }

        // time-dependent scaling factor to prevent
        // excessive growth of the numerical values of
        // alpha and beta matricies
        vec denom(len); denom.fill(-inf);

        // branch transition likelihood matricies
        mat gammau(TrellisStructure.numInputSymbols,len);
        mat gammac(TrellisStructure.numOutputSymbols,len);

        // calculate uncoded transition probability values
        {

            for (size_t i=0; i<len; i++)
            {
            
                // compute distance metric for all input combinations
                for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                {

                    mat tmp = enc.bininp.row(j) * flipud(Lu(span(i*enc.k,(i+1)*enc.k),0));
                    gammau(j,i) = 0.5 * tmp(0);

                }

            }

        }

        // calculate coded transition probability values
        {

            for (size_t i=0; i<len; i++)
            {
            
                // compute distance metric for all input combinations
                for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
                {

                    mat tmp = enc.next.binout.row(j) *
                        flipud(Lc(span(i*enc.n,(i+1)*enc.n-1),
                        0));
                    gammac(j,i) = 0.5 * tmp(0);

                }

            }

        }

        // calculate alpha values [forward recursion]
        {

            for (size_t i=0; i<len; i++)
            {

                for (size_t s=0; s<numStates; s++)
                {

                    pstate_t & pState = pStateMap[s];

                    if (pState.count == 0)
                    {
                        ak(numStates,i) = -inf;
                        continue;
                    }

                    ak(numStates,i) = gammau(pState.inpt[0],i) + gammac(pState.outp[0],i) + ak(pState.prev[0],i);

                    for (size_t j=1; j<pState.count; j++)
                    {
                        double mtr = gammau(pState.inpt[j],i) + gammac(pState.outp[j],i) + ak(pState.prev[j],i);
                        ak(s,i+1) = APPDecoderPrivate::max_(mtr,ak(s,i+1));
                    }

                    denom(i) = APPDecoderPrivate::max_(ak(s,i+1),denom(i));

                }

                // normalize alpha
                for (size_t s=0; s<numStates; s++)
                {
                    ak(s,i+1) -= denom(i);
                }

            }

        }

        // calculate beta values [backwards recursion]
        {

            for (size_t i=len; i>0; i--)
            {

                for (size_t s=0; s<numStates; s++)
                {

                    bk(s,i-1) = gammau(0,len-1) + gammac(TrellisStructure.outputs(s,0)) + bk(TrellisStructure.nextStates(s,0),len);

                    for (size_t j=1; j<TrellisStructure.numInputSymbols; j++)
                    {
                        double mtr = gammau(j,len-1) + gammac(TrellisStructure.outputs(s,j)) + bk(TrellisStructure.nextStates(s,j),len);
                        bk(s,i-1) = APPDecoderPrivate::max_(mtr,bk(s,i-1));
                    }

                }

                // normalize beta
                for (size_t s = 0; s<numStates; s++)
                {
                    bk(s,i-1) -= denom(i-1);
                }

            }

        }

        // pre-allocate unencoded bit probability
        // output vector L(u)
        lud->set_size(Lu.n_rows,Lu.n_cols);

        // update unencoded bit probability
        // log-likelihood ratio
        {

            mat & LUD = *lud;

            for (size_t i=0; i<len; i++)
            {

                for (size_t k=0; k<enc.k; k++)
                {
                    
                    float one_val = -inf;
                    float zero_val = -inf;

                    int mask = (1<<(enc.k-k-1));

                    for (size_t s=0; s<numStates; s++)
                    {

                        for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                        {

                            // calculate path metric
                            int nxt = TrellisStructure.nextStates(s,j);
                            int idx = TrellisStructure.outputs(s,j);
                            double metr = ak(s,i) + bk(nxt,i+1) + gammac(idx,i);

                            if (j&mask)     // binary one input
                                one_val = APPDecoderPrivate::max_(metr,one_val);
                            else
                                zero_val = APPDecoderPrivate::max_(metr,zero_val);

                        }

                    }

                    // update output metric
                    LUD(i*enc.k+k) = one_val - zero_val;

                }

            }

        }

        // update codeded bit probability
        // log-likelihood ratio
        if ( lcd != 0 )
        {

            // resize encoeded bit probability
            // return vector
            lcd->set_size(Lc.n_rows,Lc.n_cols);

            // convert object pointer to object reference
	        mat & LCD = *lcd;	

            for (size_t i=0; i<len; i++)
            {
                
                for (size_t n=0; n<enc.n; n++)
                {
                    
                    double one_val = -inf;
                    double zero_val = -inf;

                    for (size_t s=0; s<numStates; s++)
                    {

                        for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                        {

                            // calculate path metric
                            int nxt = TrellisStructure.nextStates(s,j);
                            int idx = TrellisStructure.outputs(s,j);
                            double metr = ak(s,i) + bk(nxt,i+1) + gammau(j,i);

                            if (idx&(1<<(enc.n-n-1)))       // binary one input
                                one_val = APPDecoderPrivate::max_(metr,one_val);
                            else                            // binary zero input
                                zero_val = APPDecoderPrivate::max_(metr,zero_val);

                        }

                    }

                    // update output metric
                    LCD(i*enc.n+n) = one_val - zero_val;

                }

            }

        }

    }

    /*!
      Construct default decoder object
    */
    APPDecoder::APPDecoder()
        : balletObject(*new APPDecoderPrivate)
    {

        // default trellis
        imat constraints = "7";
        imat codegen = "171 133;";
        TrellisStructure = poly2trellis(constraints,codegen);

        // default algorithm
        Algorithm = "Max*";

        // misc. default properties
        NumScalingBits = 3;
        TerminationMethod = "Terminated";

    }

    /*!
      Construct APPDecoder object with TrellisStructure
      set to TRELLIS
    */
    APPDecoder::APPDecoder(const Trellis &TRELLIS)
        : balletObject(*new APPDecoderPrivate)
    {
       
        // set TrellisStructure to TRELLIS 
        TrellisStructure.numInputSymbols = TRELLIS.numInputSymbols;
        TrellisStructure.numOutputSymbols = TRELLIS.numOutputSymbols;
        TrellisStructure.numStates = TRELLIS.numStates;
        TrellisStructure.nextStates = TRELLIS.nextStates;
        TrellisStructure.outputs = TRELLIS.outputs;

        // set default algorithm
        Algorithm = "Max*";

        // misc. default properties
        NumScalingBits = 3;
        TerminationMethod = "Terminated";

    }

    /*!
      Decode convolutional code using the a posteriori probability method
    */
    mat APPDecoder::decode(const mat &Lu, const mat &Lc, mat *lcd)
    {

        BALLET_D(APPDecoder);

        // lock nontunable properties and input specifications
        if (!isLocked()) { d->lock(); }

        mat y;
        d->decode(Lu,Lc,lcd,&y);

        return y;

    }

};
