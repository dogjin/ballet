#include "lva_p.h"

namespace ballet
{

    ListViterbiResultPrivate::ListViterbiResultPrivate()
        : balletObjectPrivate()
          , objectData(0)
          , serialObjectData(0)
          , num_paths(0)
          , current(0)
    {
    }

    ListViterbiResultPrivate::~ListViterbiResultPrivate()
    {

        if (objectData)
            delete objectData;
        if (serialObjectData)
            delete serialObjectData;
    }

    ListViterbiDecoderPrivate::ListViterbiDecoderPrivate()
        : balletObjectPrivate()
          , extra(0)
    {
    }

    ListViterbiDecoderPrivate::~ListViterbiDecoderPrivate()
    {
        if (extra)
            deleteExtra();
    }

    void ListViterbiDecoderPrivate::createExtra()
    {

        if (!extra)
        {
            extra = new VDExtra;
        }

    }

    void ListViterbiDecoderPrivate::deleteExtra()
    {
        if (extra)
        {
            delete extra;
        }
    }

    void ListViterbiDecoderPrivate::lock()
    {
        BALLET_Q(ListViterbiDecoder);

        if (!extra)
            createExtra();
        
        // lock TrellisStructure property
        TrellisStructure = q->TrellisStructure;

        // lock InputFormat property
        InputFormat = q->InputFormat;

        // lock SoftInputWordLength
        SoftInputWordLength = q->SoftInputWordLength;

        // lock TerminationMethod property
        TerminationMethod = q->TerminationMethod;

        // lock OutputFormat property
        OutputFormat = q->OutputFormat;

        // lock Algorithm property
        Algorithm = q->Algorithm;

        // lock L property
        L = q->L;

        // number of bits per symbol
        double kdbl = std::log((double)TrellisStructure.numInputSymbols) / std::log(2.0);
        extra->k = static_cast<size_t>(kdbl);
        double ndbl = std::log((double)TrellisStructure.numOutputSymbols) / std::log(2.0);
        extra->n = static_cast<size_t>(ndbl);

        // number of tail bits
        double numTails_ = log((double)TrellisStructure.numStates)/log(2.0);
        numTails = static_cast<int>(numTails_+0.5);

        // calculate bit output values
        extra->binout.set_size(TrellisStructure.numOutputSymbols * extra->n);
        if (InputFormat.compare("Unquantized") == 0)
        {
            for (size_t i=0; i<TrellisStructure.numOutputSymbols; i++)
            {
                int out_symbol = i;
                for (size_t j=0; j<extra->n; j++)
                {
                    extra->binout(i*extra->n+j) = -2.0*static_cast<double>(out_symbol&1) + 1.0;
                    out_symbol = out_symbol >> 1;
                }
            }
        }
        else if (InputFormat.compare("Hard") == 0)
        {
            for (size_t i=0; i<TrellisStructure.numOutputSymbols; i++)
            {
                int out_symbol = i;
                for (size_t j=0; j<extra->n; j++)
                {
                    extra->binout(i*extra->n+j) = static_cast<double>(out_symbol&1);
                    out_symbol = out_symbol >> 1;
                }
            }
        }
        else if (InputFormat.compare("Soft") == 0)
        {

            extra->soft_zero = 0.0;
            extra->soft_one = static_cast<double>((2<<SoftInputWordLength)-1);

            for (size_t i=0; i<TrellisStructure.numOutputSymbols; i++)
            {
                int out_symbol = i;
                for (size_t j=0; j<extra->n; j++)
                {
                    extra->binout(i*extra->n+j) = extra->soft_one * static_cast<double>(out_symbol&1);
                    out_symbol = out_symbol >> 1;
                }
            }

        }

        // set property_locked flag
        property_locked = true;
        
    }

    ListViterbiResult ListViterbiDecoderPrivate::decode(const splib::fvec &x)
    {

        if (Algorithm.compare("Parallel") == 0)
        {
            return decode_parallel(x);
        }
        else if (Algorithm.compare("Serial") == 0)
        {
            return decode_serial(x);
        }

        return ListViterbiResult();

    }

    ListViterbiResult ListViterbiDecoderPrivate::decode_parallel(const splib::fvec &x)
    {

        ListViterbiResult decodeResult;

        // define infinity value
        double inf = numeric_limits<double>::max();

        // number of trellis states
        size_t numStates = TrellisStructure.numStates;

        // number number of trellis sections
        size_t LEN = x.length() / extra->n;

        // general decoder object data
        ListViterbiPrivate::ObjectData * objectData = new ListViterbiPrivate::ObjectData;
        objectData->TrellisStructure = TrellisStructure;
        objectData->M = LEN;
        objectData->L = L;
        objectData->path_matrix.resize(L);
        objectData->k = extra->k;
        objectData->numTails = numTails;
        objectData->OutputFormat = OutputFormat;

        // [branch metric unit]
        if (InputFormat.compare("Unquantized") == 0)
        {
            objectData->branchmetric = branch_metric_unq(x);
        }
        else if (InputFormat.compare("Hard") == 0)
        {
            objectData->branchmetric = branch_metric_har(x);
        }
        else if (InputFormat.compare("Soft") == 0)
        {
            objectData->branchmetric = branch_metric_sof(x);
        }

        // pointer to trellis members
        const int * NEXT = TrellisStructure.nextStates.begin();
        const int * OUT = TrellisStructure.codeOutputs.begin();

        // path metric unit memory allocation
        splib::vec phi((LEN+1)*numStates*L);
        splib::ivec xi(LEN*numStates*L);
        splib::ivec lmbda(LEN*numStates*L);

        // reference to trellis computation variables
        splib::vec& branchmetric = objectData->branchmetric;

        // [path metric unit]
        {

            // [initialization]
            phi = -inf;
            phi(0) = 0;

            for (size_t i=0; i<LEN; i++)
            {

                // initialize trellis transition pointers
                double *GAMMAC = branchmetric.begin() + i*TrellisStructure.numOutputSymbols;
                double * PHI = phi.begin() + i*numStates*L;
                int * LAMBDA = lmbda.begin() + i*numStates*L;
                int * XI = xi.begin() + i*numStates*L;
                
                for (size_t s=0; s<numStates; s++)
                {
                    for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                    {

                        // calculate path metric
                        int nxt = NEXT[j*numStates+s];
                        int idx = OUT[j*numStates+s];
                        double gammac = GAMMAC[idx];

                        for (size_t l=0; l<L; l++)
                        {

                            double mtr = gammac + PHI[L*s+l];

                            // update and store path metric
                            for (size_t k=0; k<L; k++)
                            {

                                if (mtr > PHI[L*(numStates+nxt)+k])
                                {

                                    // shift in place
                                    for (size_t m=L-1; m!=k; m--)
                                    {
                                        PHI[L*(numStates+nxt)+m] = PHI[L*(numStates+nxt)+m-1];
                                        LAMBDA[L*nxt+m] = LAMBDA[L*nxt+m-1];
                                        XI[L*nxt+m] = XI[L*nxt+m-1];
                                    }
                                        
                                    PHI[L*(numStates+nxt)+k] = mtr;
                                    LAMBDA[L*nxt+k] = j;
                                    XI[L*nxt+k] = L*s + l;
                                    break;
    
                                }

                            }

                        }

                    }

                }

            }

        }

        // [traceback]
        {

            for (size_t l=0; l<L; l++)
            {

                // [termination state]
                int curstate = l;

                // initialize return vector
                std::vector<ListViterbiPrivate::path_t> mypath(LEN+1);

                // trellis termination for l:th best path
                mypath[LEN].state = 0;
                mypath[LEN].phi = phi(LEN*numStates*L+l);
                mypath[LEN].chi = 0.0;
                mypath[LEN].lmbda = -1;

                for (size_t i=LEN; i>0; i--)
                {

                    // initialize trellis section pointers
                    double * PHI = phi.begin() + (i-1)*numStates*L;
                    int * LAMBDA = lmbda.begin() + (i-1)*numStates*L;
                    int * XI = xi.begin() + (i-1)*numStates*L;

                    // most likely input symbol
                    mypath[i-1].state = XI[curstate]/L;
                    mypath[i-1].phi = PHI[XI[curstate]];
                    mypath[i-1].chi = mypath[LEN].phi - PHI[XI[curstate]];
                    mypath[i-1].lmbda = LAMBDA[curstate];

                    // traceback to previous trellis stage
                    curstate = XI[curstate];

                }

                objectData->path_matrix[l] = mypath;

            }

        }

        // finalize result parameters
        ListViterbiResultPrivate * d = decodeResult.d_func();
        d->num_paths = L;
        d->objectData = objectData;

        //return decode result
        return decodeResult;

    }

    ListViterbiResult ListViterbiDecoderPrivate::decode_serial(const splib::fvec &x)
    {

        ListViterbiResult decodeResult;

        // define infinity value
        double inf = numeric_limits<double>::max();

        // number of trellis states
        size_t numStates = TrellisStructure.numStates;

        // number number of trellis sections
        size_t LEN = x.length() / extra->n;

        // general decoder object data
        ListViterbiPrivate::ObjectData * objectData = new ListViterbiPrivate::ObjectData;
        objectData->TrellisStructure = TrellisStructure;
        objectData->M = LEN;
        objectData->L = L;
        objectData->path_matrix.resize(L);
        objectData->k = extra->k;
        objectData->numTails = numTails;
        objectData->OutputFormat = OutputFormat;

        // serial decoder object data
        ListViterbiPrivate::SerialObjectData * serialObjectData = new ListViterbiPrivate::SerialObjectData;
        serialObjectData->merge_count.set_size(LEN);
        serialObjectData->merge_count.ones();
        serialObjectData->merge_matrix.set_size(LEN,L);
        serialObjectData->merge_matrix.zeros();

        // [branch metric unit]
        if (InputFormat.compare("Unquantized") == 0)
        {
            objectData->branchmetric = branch_metric_unq(x);
        }
        else if (InputFormat.compare("Hard") == 0)
        {
            objectData->branchmetric = branch_metric_har(x);
        }
        else if (InputFormat.compare("Soft") == 0)
        {
            objectData->branchmetric = branch_metric_sof(x);
        }

        // pointer to trellis members
        const int * NEXT = TrellisStructure.nextStates.begin();
        const int * OUT = TrellisStructure.codeOutputs.begin();

        // path metric unit memory allocation
        serialObjectData->phi.set_size((LEN+1)*numStates);
        serialObjectData->chi.set_size((LEN+1)*numStates);
        serialObjectData->xi.set_size(LEN*numStates);
        serialObjectData->lmbda.set_size(LEN*numStates);

        // reference to trellis computation variables
        splib::vec& branchmetric = objectData->branchmetric;
        splib::vec& phi = serialObjectData->phi;
        splib::ivec& xi = serialObjectData->xi;
        splib::ivec& lmbda = serialObjectData->lmbda;

        // [path metric unit]
        {

            // [initialization]
            phi = -inf;
            phi(0) = 0;
            lmbda.zeros();

            for (size_t i=0; i<LEN; i++)
            {

                // initialize trellis transition pointers
                double *GAMMAC = branchmetric.begin() + i*TrellisStructure.numOutputSymbols;
                double * PHI = phi.begin() + i*numStates;
                int * LAMBDA = lmbda.begin() + i*numStates;
                int * XI = xi.begin() + i*numStates;
                
                for (size_t s=0; s<numStates; s++)
                {
                    for (size_t j=0; j<TrellisStructure.numInputSymbols; j++)
                    {

                        // calculate path metric
                        int nxt = NEXT[j*numStates+s];
                        int idx = OUT[j*numStates+s];
                        double mtr = GAMMAC[idx] + PHI[s];

                        // update and store path metric
                        if (mtr > PHI[numStates+nxt])
                        {
                            PHI[numStates+nxt] = mtr;
                            LAMBDA[nxt] = j;
                            XI[nxt] = s;
                        }

                    }

                }

            }

        }

        // [traceback]
        {

            // allocate array for best path 
            // through trellis
            std::vector<ListViterbiPrivate::path_t> best_path(LEN+1);

            // [termination state]
            int curstate = 0;

            // trellis termination for best path
            // through trellis
            best_path[LEN].state = curstate;
            best_path[LEN].phi = phi(LEN*numStates);
            best_path[LEN].chi = 0.0;
            best_path[LEN].lmbda = -1;

            for (size_t i=LEN; i>0; i--)
            {

                // initialize trellis section pointers
                double * PHI = phi.begin() + (i-1)*numStates;
                int * LAMBDA = lmbda.begin() + (i-1)*numStates;
                int * XI = xi.begin() + (i-1)*numStates;

                // most likely path
                best_path[i-1].state = XI[curstate];
                best_path[i-1].phi = PHI[XI[curstate]];
                best_path[i-1].chi = best_path[LEN].phi-best_path[i-1].phi;
                best_path[i-1].lmbda = LAMBDA[curstate];

                // traceback to previous trellis stage
                curstate = XI[curstate];

            }

            // add globally best path to path_matrix
            objectData->path_matrix[0] = best_path;

        }

        // finalize result parameters
        ListViterbiResultPrivate * d = decodeResult.d_func();
        d->num_paths = 1;
        d->objectData = objectData;
        d->serialObjectData = serialObjectData;

        return decodeResult;

    }

    splib::ivec ListViterbiResult::operator()(const int n)
    {
        BALLET_D(ListViterbiResult);
        return (*d)(n);
    }

    splib::ivec ListViterbiResult::next()
    {

        BALLET_D(ListViterbiResult);
        return d->next();

    }

    typedef struct
    {
        int final_global_merge;
        int local_rank;
        int n;
        int final_local_merge;
        int state;
        int input;
    } nextpath_t;

    splib::ivec ListViterbiResultPrivate::next()
    {

        if ( num_paths<=current )
        {

            // define infinity value
            double inf = numeric_limits<double>::max();

            // number of trellis states
            size_t numStates = objectData->TrellisStructure.numStates;

            // number number of trellis sections
            size_t LEN = objectData->M;

            // get reference to best path through trellis
            std::vector<ListViterbiPrivate::path_t> & best_path = objectData->path_matrix[0];

            // pointer to trellis members
            const int * NEXT = objectData->TrellisStructure.nextStates.begin();
            const int * OUT = objectData->TrellisStructure.codeOutputs.begin();

            while ( num_paths<=current )
            {

                // initialize k:th best candidate
                double nxtscore = -inf;
                nextpath_t nxtpath = {0,0,0,0,0,0};

                // get score of globally n:th best path through the trellis
                double cscore = objectData->path_matrix[num_paths-1].back().phi;

                for (size_t t=0; t<LEN; t++)
                {

                    // initialize trellis transition pointer
                    double * PHI = serialObjectData->phi.begin() + t*numStates;
                    double * GAMMAC = objectData->branchmetric.begin() + t*objectData->TrellisStructure.numOutputSymbols;

                    // get number of paths that finally
                    // merge with the globally best path
                    // at time t
                    int c = serialObjectData->merge_count[t];

                    // state occupied by best path at time t
                    int cur = best_path[t+1].state;

                    if (c == 1)
                    {

                        // cost remaining by following best
                        // path for the remainder of the trellis
                        double remaining_cost = best_path[t+1].chi;

                        // previous state by following the best path    
                        int prev = best_path[t].state;

                        for (size_t j=0; j<numStates; j++)
                        {
                            for (size_t nput=0; nput<objectData->TrellisStructure.numInputSymbols; nput++)
                            {
                                
                                // avoid tracing path of best_path
                                if ( j == prev ) { continue; }
                                
                                // state transition j -> i
                                int nxt = NEXT[nput*numStates+j];
                                if ( nxt != cur ) { continue; }

                                // transition output, transition metric
                                int idx = OUT[nput*numStates+j];
                                double gammac = GAMMAC[idx];

                                // update and store path metric
                                double mtr = gammac + PHI[j];
                                if ( !(mtr+remaining_cost < cscore) ) { continue; }
                                if ( mtr+remaining_cost > nxtscore )
                                {
                                    nxtscore = mtr+remaining_cost;
                                    nxtpath.final_global_merge = t;
                                    nxtpath.local_rank = serialObjectData->merge_count[t]+1;
                                    nxtpath.n = 0;
                                    nxtpath.final_local_merge = t;
                                    nxtpath.state = j;
                                    nxtpath.input = nput;
                                }
                            
                            }

                        }

                    }
                    else
                    {

                        for (size_t k=0; k<c; k++)
                        {

                            // state occupied by best path at time t
                            int cur = best_path[t+1].state;

                            // get pointer to c:th best path that
                            // finally merges at time t
                            int n = serialObjectData->merge_matrix(t,c-k-1);
                            std::vector<ListViterbiPrivate::path_t> & cpath = objectData->path_matrix[n];

                            for (size_t rtau=0; rtau<=t; rtau++)
                            {

                                size_t tau = t-rtau;

                                // trellis transition pointers
                                PHI = serialObjectData->phi.begin() + tau*numStates;
                                GAMMAC = objectData->branchmetric.begin() + tau*objectData->TrellisStructure.numOutputSymbols;

                                // cost remaining by following c:th best
                                // path for the remainder of the trellis
                                double remaining_cost = cpath[tau+1].chi;
                                
                                // previous state by following c:th best path
                                int prev = cpath[tau].state;
                    
                                for (size_t j=0; j<numStates; j++)
                                {
                                    for (size_t nput=0; nput<objectData->TrellisStructure.numInputSymbols; nput++)
                                    {
                
                                        // avoid tracing path of c
                                        if ( j == prev ) { continue; }
                
                                        // state transition j -> i
                                        int nxt = NEXT[nput*numStates+j];
                                        if ( nxt != cur ) { continue; }
                
                                        // transition output, transition metric
                                        int idx = OUT[nput*numStates+j];
                                        double gammac = GAMMAC[idx];
                
                                        // update and store path metric
                                        double mtr = gammac + PHI[j];
                                        if ( !(mtr+remaining_cost < cscore) ) { continue; }
                                        if ( mtr+remaining_cost > nxtscore )
                                        {
                                            nxtscore = mtr+remaining_cost;
                                            nxtpath.final_global_merge = t;
                                            nxtpath.local_rank = serialObjectData->merge_count[t]+1;
                                            nxtpath.n = n;
                                            nxtpath.final_local_merge = tau;
                                            nxtpath.state = j;
                                            nxtpath.input = nput;
                                        }

                                    }

                                }

                                // follow c:th best path to 
                                // beginning of the trellis
                                cur = prev;

                            }

                        }

                    }

                }

                {

                    // time instant of final global merge
                    size_t t = nxtpath.final_global_merge;

                    // get number of paths that finally
                    // merge with the globally best path
                    // at time t
                    int c = serialObjectData->merge_count[t];

                    // get pointer to c:th best path that
                    // finally merges at time t
                    std::vector<ListViterbiPrivate::path_t> & next_best = objectData->path_matrix[nxtpath.n];
                    std::vector<ListViterbiPrivate::path_t> mypath(LEN+1);

                    // [pre-merge]
                    int s = nxtpath.state;
                    int flm = nxtpath.final_local_merge;
                    for (size_t rtau=0; rtau<=flm; rtau++)
                    {

                        size_t tau = flm-rtau;

                        // trellis transition pointers
                        double * PHI = serialObjectData->phi.begin() + tau*numStates;

                        // traceback best path from state s
                        // to beginning of trellis
                        mypath[tau].state = s;
                        mypath[tau].phi = PHI[s];

                        // get previous state
                        if (tau != 0)
                        {
                            mypath[tau-1].lmbda = serialObjectData->lmbda((tau-1)*numStates+s);
                            s = serialObjectData->xi((tau-1)*numStates+s);
                        }

                    }

                    // [merge]
                    int idx = OUT[nxtpath.input*numStates+nxtpath.state];
                    mypath[flm+1].state = next_best[flm+1].state;
                    mypath[flm+1].phi = mypath[flm].phi + objectData->branchmetric(flm*objectData->TrellisStructure.numOutputSymbols+idx);
                    mypath[flm].lmbda = nxtpath.input;
                    double diff = mypath[flm+1].phi - next_best[flm+1].phi;

                    // [post-merge]
                    for (size_t tau=flm+2; tau<LEN+1; tau++)
                    {

                        // remainder same as globally next better path
                        mypath[tau].state = next_best[tau].state;
                        mypath[tau].phi = next_best[tau].phi + diff;
                        mypath[tau-1].lmbda = next_best[tau-1].lmbda;

                    }

                    // compute remainder
                    for (size_t tau=0; tau<LEN+1; tau++)
                    {
                        mypath[tau].chi = mypath[LEN].phi - mypath[tau].phi;
                    }

                    // update merge count
                    serialObjectData->merge_count(t) = serialObjectData->merge_count(t) + 1;

                    // add path to matrix
                    objectData->path_matrix[num_paths] = mypath;
                    serialObjectData->merge_matrix(t,c) = num_paths;
                    num_paths++;

                }

            }

        }

        // allocate return array
        splib::ivec symbols(objectData->M);

        // current next-best vector
        std::vector<ListViterbiPrivate::path_t> & mypath = objectData->path_matrix[current];
       
        // convert path to decoded sequence
        for (size_t n=0; n<objectData->M; n++)
        {
            symbols(n) = mypath[n].lmbda;
        }

        // translate decoded sequence to proper output format
        splib::ivec y;
        if (objectData->OutputFormat.compare("Binary") == 0)
        {

            // allocate output array
            y.set_size((objectData->M-objectData->numTails) * objectData->k);
            for (size_t n=0; n<objectData->M-objectData->numTails; n++)
            {

                // initialize pointer [fast access]
                int * Y = y.begin() + n*objectData->k;

                // decoded symbol value
                int input_symbol = symbols(n);
   
                // convert symbol to binary vector 
                for (size_t k=0; k<objectData->k; k++)
                {
                    Y[objectData->k-k-1] = static_cast<int>(input_symbol&1);
                    input_symbol = input_symbol >> 1;
                }
                
            }

        }
        else if (objectData->OutputFormat.compare("Integer") == 0)
        {
    
            // copy decoded symbols to output array
            y = symbols(0,objectData->M-objectData->numTails-1);

        }

        // update read pointer
        current++;

        return y;

    }

    splib::ivec ListViterbiResultPrivate::operator()(const int n)
    {

        while (num_paths <= n) { next(); }

        // allocate return array
        splib::ivec symbols(objectData->M);

        // n:th best vector
        std::vector<ListViterbiPrivate::path_t> & mypath = objectData->path_matrix[n];
       
        // convert path to decoded sequence
        for (size_t i=0; i<objectData->M; i++)
        {
            symbols(i) = mypath[i].lmbda;
        }

        // translate decoded sequence to proper output format
        splib::ivec y;
        if (objectData->OutputFormat.compare("Binary") == 0)
        {

            // allocate output array
            y.set_size((objectData->M-objectData->numTails) * objectData->k);
            for (size_t i=0; i<objectData->M-objectData->numTails; i++)
            {

                // initialize pointer [fast access]
                int * Y = y.begin() + i*objectData->k;

                // decoded symbol value
                int input_symbol = symbols(i);
   
                // convert symbol to binary vector 
                for (size_t k=0; k<objectData->k; k++)
                {
                    Y[objectData->k-k-1] = static_cast<int>(input_symbol&1);
                    input_symbol = input_symbol >> 1;
                }
                
            }

        }
        else if (objectData->OutputFormat.compare("Integer") == 0)
        {
    
            // copy decoded symbols to output array
            y = symbols(0,objectData->M-objectData->numTails-1);

        }

        // update current pointer
        current = n + 1;

        return y;

    }

    splib::vec ListViterbiDecoderPrivate::branch_metric_unq(const splib::fvec &x)
    {

        size_t LEN = x.length() / extra->n;

        // pointer to internal data structures
        const double * BINOUT = extra->binout.begin();

        // allocate branch metric array
        splib::vec y(LEN * TrellisStructure.numOutputSymbols);

        // [branch metric]
        for (size_t i=0; i<LEN; i++)
        {

            // get output symbol i
            const float * buf = x.begin() + i*extra->n;

            // compute distance metric for all input combinations
            for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
            {

                double tmp = 0.0;
                for (size_t n=0; n<extra->n; n++)
                {
                    tmp = tmp + static_cast<double>(buf[extra->n-n-1]) * BINOUT[extra->n*j + n];
                }

                y(i*TrellisStructure.numOutputSymbols+j) = tmp;

            }

        }

        return y;

    }

    splib::vec ListViterbiDecoderPrivate::branch_metric_har(const splib::fvec &x)
    {

        size_t LEN = x.length() / extra->n;

        // pointer to internal data structures
        const double * BINOUT = extra->binout.begin();

        // allocate branch metric array
        splib::vec y(LEN * TrellisStructure.numOutputSymbols);

        // [branch metric]
        for (size_t i=0; i<LEN; i++)
        {

            // get output symbol i
            const float * buf = x.begin() + i*extra->n;

            // compute distance metric for all input combinations
            for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
            {

                double tmp = 0.0;
                for (size_t n=0; n<extra->n; n++)
                {
                    tmp = tmp + static_cast<double>(~(static_cast<int>(buf[extra->n-n-1])^static_cast<int>(BINOUT[extra->n*j + n])));
                }

                y(i*TrellisStructure.numOutputSymbols+j) = tmp;

            }

        }

        return y;

    }

    splib::vec ListViterbiDecoderPrivate::branch_metric_sof(const splib::fvec &x)
    {

        size_t LEN = x.length() / extra->n;

        // pointer to internal data structures
        const double * BINOUT = extra->binout.begin();

        // allocate branch metric array
        splib::vec y(LEN * TrellisStructure.numOutputSymbols);

        // [branch metric]
        for (size_t i=0; i<LEN; i++)
        {

            // get output symbol i
            const float * buf = x.begin() + i*extra->n;

            // compute distance metric for all input combinations
            for (size_t j=0; j<TrellisStructure.numOutputSymbols; j++)
            {

                double tmp = 0.0;
                for (size_t n=0; n<extra->n; n++)
                {
                    tmp = tmp + extra->soft_one - std::abs(static_cast<double>(buf[extra->n-n-1])-BINOUT[extra->n*j + n]);
                }

                y(i*TrellisStructure.numOutputSymbols+j) = tmp;

            }

        }

        return y;

    }

    ListViterbiResult::ListViterbiResult()
        : balletObject(*new ListViterbiResultPrivate)
    {
    }

    ListViterbiResult::ListViterbiResult(const ListViterbiResult &other)
        : balletObject(*new ListViterbiResultPrivate)
    {

        // pointer to private class
        BALLET_D(ListViterbiResult);
        
        // pointer to other private class
        const ListViterbiResultPrivate * other_d = other.d_func();

        // copy object data
        ListViterbiPrivate::ObjectData * objectData = new ListViterbiPrivate::ObjectData;
        objectData->TrellisStructure = other_d->objectData->TrellisStructure;
        objectData->M = other_d->objectData->M;
        objectData->L = other_d->objectData->L;
        objectData->path_matrix.resize(objectData->L);
        objectData->k = other_d->objectData->k;
        objectData->OutputFormat = other_d->objectData->OutputFormat;
        objectData->branchmetric = other_d->objectData->branchmetric;

        // copy serial object data [if applicable]
        ListViterbiPrivate::SerialObjectData * serialObjectData = NULL;
        if ( other_d->serialObjectData )
        {

            serialObjectData = new ListViterbiPrivate::SerialObjectData;
            serialObjectData->phi = other_d->serialObjectData->phi;
            serialObjectData->chi = other_d->serialObjectData->chi;
            serialObjectData->xi = other_d->serialObjectData->xi;
            serialObjectData->lmbda = other_d->serialObjectData->lmbda;
            serialObjectData->merge_count.set_size(objectData->M);
            serialObjectData->merge_count.ones();
            serialObjectData->merge_matrix.set_size(objectData->M,objectData->L);
            serialObjectData->merge_matrix.zeros();
            
        }

        // copy number of paths in path_matrix 
        d->num_paths = other_d->num_paths;

        // copy paths from path_matrix
        for (size_t n=0; n<d->num_paths; n++)
            objectData->path_matrix[n] = other_d->objectData->path_matrix[n];

        // finalize object
        d->objectData = objectData;
        d->serialObjectData = serialObjectData;
        
    }

    ListViterbiDecoder::ListViterbiDecoder()
        : balletObject(*new ListViterbiDecoderPrivate)
    {

        // default trellis
        splib::ivec constraints = "7";
        splib::imat codegen = "177 133;";
        TrellisStructure = splib::poly2trellis(constraints,codegen);

        // default InputFormat
        InputFormat = "Unquantized";

        // default OutputFormat
        OutputFormat = "Integer";

        // default Algorithm
        Algorithm = "Parallel";

        // default L value
        L = 4;

    }

    ListViterbiResult & ListViterbiResult::operator=(const ListViterbiResult &other)
    {

        if (this == &other) { return *this; }

        // pointer to private class
        BALLET_D(ListViterbiResult);
        
        // pointer to other private class
        const ListViterbiResultPrivate * other_d = other.d_func();

        // copy object data [if applicable]
        if ( other_d->objectData == NULL )
        {
            if ( d->objectData != NULL )
            {
                delete d->objectData;
                d->objectData = NULL;
            }
        }
        else
        {
            if ( d->objectData == NULL )
                d->objectData = new ListViterbiPrivate::ObjectData;
            d->objectData->TrellisStructure = other_d->objectData->TrellisStructure;
            d->objectData->M = other_d->objectData->M;
            d->objectData->L = other_d->objectData->L;
            d->objectData->path_matrix.resize(d->objectData->L);
            d->objectData->k = other_d->objectData->k;
            d->objectData->OutputFormat = other_d->objectData->OutputFormat;
            d->objectData->branchmetric = other_d->objectData->branchmetric;
        }

        // copy serial object data [if applicable]
        if ( other_d->serialObjectData )
        {

            if ( d->serialObjectData == NULL )
                d->serialObjectData = new ListViterbiPrivate::SerialObjectData;

            d->serialObjectData->phi = other_d->serialObjectData->phi;
            d->serialObjectData->chi = other_d->serialObjectData->chi;
            d->serialObjectData->xi = other_d->serialObjectData->xi;
            d->serialObjectData->lmbda = other_d->serialObjectData->lmbda;
            d->serialObjectData->merge_count.set_size(d->objectData->M);
            d->serialObjectData->merge_count.ones();
            d->serialObjectData->merge_matrix.set_size(d->objectData->M,d->objectData->L);
            d->serialObjectData->merge_matrix.zeros();
            
        }

        // copy number of paths in path_matrix 
        d->num_paths = other_d->num_paths;

        // copy paths from path_matrix
        for (size_t n=0; n<d->num_paths; n++)
            d->objectData->path_matrix[n] = other_d->objectData->path_matrix[n];

        return *this;

    }


    ListViterbiResult ListViterbiDecoder::decode(const splib::fvec &x)
    {

        BALLET_D(ListViterbiDecoder);

        // lock decoder object
        if (!isLocked()) { d->lock(); }

        return d->decode(x);

    }


};
