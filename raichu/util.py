import numpy as np
from collections import Counter
from joblib import Parallel, delayed
import cooler, h5py, numba, logging
from scipy.optimize import dual_annealing
from numba import njit, double, int32
from collections import defaultdict
        
@njit(double(double[:], numba.types.Tuple((double[:], int32[:,:], double[:]))), parallel=True, fastmath=True, nogil=True, cache=True)
def eval_func(weights, *args):
    
    data, coords, Earr = args
    weights = np.power(10, weights)
    W = np.multiply(weights[coords[:,0]], weights[coords[:,1]])
    obj_arr = np.abs(data - Earr * W)
    obj = obj_arr.sum()

    return obj

def optimize_by_dual_annealing(ini_weights, data, coords, Earr, lb, ub, maxiter):

    # data and coords are returned by "extract_valid_pixels"
    lw = [lb] * ini_weights.size
    up = [ub] * ini_weights.size
    
    # optimize the weights using dual annealing
    ret = dual_annealing(
        eval_func,
        bounds=list(zip(lw, up)),
        args=(data, coords, Earr),
        maxiter=maxiter,
        seed=42,
        x0=ini_weights,
        no_local_search=False,
        minimizer_kwargs={'method':'L-BFGS-B'}
    )

    weights = ret['x']
    
    return weights, ret

def calculate_expected_core(clr, c, included_bins, max_dis):

    M = clr.matrix(balance=False, sparse=True).fetch(c).tocsr()
    marg = np.array(M.sum(axis=0)).ravel()
    valid_cols = marg > 0
    if not included_bins is None:
        tmp = np.zeros(valid_cols.size, dtype=bool)
        tmp[included_bins] = True
        valid_cols = valid_cols & tmp
    n = M.shape[0]

    expected = {}
    maxdis = min(n-1, max_dis)
    for i in range(maxdis+1):
        if i == 0:
            valid = valid_cols
        else:
            valid = valid_cols[:-i] * valid_cols[i:]
            
        diag = M.diagonal(i)[valid]
        if diag.size > 0:
            expected[i] = [diag.sum(), diag.size]
    
    return expected

def calculate_expected(clr, chroms, included_bins, max_dis, nproc=1, N=400, dynamic_window_size=10):
    
    binsize = clr.binsize
    max_dis = max_dis // binsize
    queue = []
    for c in chroms:
        if included_bins is None:
            queue.append((clr, c, None, max_dis))
        else:
            if c in included_bins:
                queue.append((clr, c, included_bins[c], max_dis))
    
    results = Parallel(n_jobs=nproc)(delayed(calculate_expected_core)(*i) for i in queue)
    diag_sums = []
    pixel_nums = []
    for i in range(max_dis+1):
        nume = 0
        denom = 0
        for extract in results:
            if i in extract:
                nume += extract[i][0]
                denom += extract[i][1]
        diag_sums.append(nume)
        pixel_nums.append(denom)
    
    Ed = {}
    for i in range(max_dis+1):
        for w in range(dynamic_window_size+1):
            tmp_sums = diag_sums[max(i-w,0):i+w+1]
            tmp_nums = pixel_nums[max(i-w,0):i+w+1]
            n_count = sum(tmp_sums)
            n_pixel = sum(tmp_nums)
            if n_count > N:
                Ed[i] = n_count / n_pixel
                break
    
    return Ed

def initialize_weights(M, included_bins):

    # input M: sparse matrix in CSR
    tmp = np.array(M.sum(axis=0)).ravel()
    if included_bins is None:
        marg = tmp
    else:
        marg = np.zeros_like(tmp)
        marg[included_bins] = tmp[included_bins]

    logNzMarg = np.log(marg[marg>0])
    med_logNzMarg = np.median(logNzMarg)
    dev_logNzMarg = cooler.balance.mad(logNzMarg)
    # MAD-max filter similar to cooler.balance
    cutoff = np.exp(med_logNzMarg - 5 * dev_logNzMarg)
    marg[marg<cutoff] = 0

    scale = marg[marg>0].mean()
    weights = np.sqrt(marg / scale)
    weights = weights.astype(np.float64)
    valid_cols = weights > 0

    return weights, valid_cols

def extract_valid_pixels(M, Ed, valid_cols, start_diag=0, top_per=99, bottom_per=1):

    # input M: sparse matrix in CSR
    # input Ed: the average contact frequency at each genomic distance
    # top_per: percentile of the upper bound O/E value
    # bottom_per: percentile of the lower bound O/E value
    x_arr = np.array([], dtype=np.int32)
    y_arr = np.array([], dtype=np.int32)
    data = np.array([], dtype=np.float64)
    OE = np.array([], dtype=np.float64)
    maxd = max(Ed)
    idx = np.arange(M.shape[0]).astype(np.int32)
    # starting from distance 0 will make the result contain more stripes
    for i in range(start_diag, maxd+1):
        diag = M.diagonal(i)
        if diag.size > 0:
            if i == 0:
                xi = idx
                yi = idx
                valid = valid_cols
            else:
                xi = idx[:-i]
                yi = idx[i:]
                valid = valid_cols[:-i] * valid_cols[i:]
            oe = diag / Ed[i]
            mask = (diag > 0) & valid
            x_arr = np.r_[x_arr, xi[mask]]
            y_arr = np.r_[y_arr, yi[mask]]
            OE = np.r_[OE, oe[mask]]
            data = np.r_[data, diag[mask]]
    
    top_v = np.percentile(OE, top_per)
    bottom_v = np.percentile(OE, bottom_per)
    mask = (OE > bottom_v) & (OE < top_v)
    x_arr = x_arr[mask]
    y_arr = y_arr[mask]
    data = data[mask]
    coords = np.r_['1,2,0', x_arr, y_arr]

    return coords, data

def pipeline(clr, chrom, Ed, ws_bin, included_bins, ndiag, lb, ub, maxiter, min_nnz, n_threads):

    # clr: Cooler
    # chrom: chromosome name
    # Ed: the average contact frequency at each genomic distance
    # ws_bin: window size in the unit of pixels
    # maxiter: the maximum number of global search iterations for dual annealing
    numba.set_num_threads(n_threads)
    logger = logging.getLogger()
    M = clr.matrix(balance=False, sparse=True).fetch(chrom).tocsr()
    indices = clr.bins().fetch(chrom).index.values
    ini_weights, valid_cols = initialize_weights(M, included_bins)
    ini_weights[ini_weights>0] = np.log10(ini_weights[ini_weights>0])
    coords, data = extract_valid_pixels(M, Ed, valid_cols, start_diag=ndiag)
    # optimize the weights in sliding windows
    ws_bp = ws_bin * clr.binsize # window size in the unit of base pairs
    if included_bins is None:
        queue = split_chromosome(clr, chrom, ws_bp)
    else:
        queue = extract_valid_regions(clr, included_bins, ws_bp)

    collect = {}
    for s, e in queue:
        #print('current region: {0}:{1}-{2}'.format(chrom, s*clr.binsize, e*clr.binsize))
        #time_start = time.time()
        ini_weights_ = ini_weights[s:e]
        mask = (coords[:,0] >= s) & (coords[:,1] < e)
        rl = e - s
        if (valid_cols[s:e].sum() / rl < 0.1) and (included_bins is None):
            # skip regions with very few data points
            weights_ = ini_weights_.copy() * np.nan
            collect[(s, e)] = [weights_, None]
            continue
        coords_ = coords[mask] - s
        data_ = data[mask]
        Earr_ = np.array([Ed[i[1]-i[0]] for i in coords_], dtype=np.float64)
        #obj_ini = eval_func(ini_weights_, data_, coords_, Earr_)
        weights_, ret = optimize_by_dual_annealing(ini_weights_, data_, coords_, Earr_,
                                                   lb, ub, maxiter)
        collect[(s, e)] = [weights_, ret]
        #elapse = time.time() - time_start
        #print('{0}s elapsed, eval func {1} -> {2}'.format(elapse, obj_ini, ret.fun))
    
    weights = combine_weights(collect, ini_weights)
    # count the number of data points for each bin
    counts = Counter(coords.ravel())
    # remove bins with <min_nnz data points
    for i in range(len(weights)):
        if counts[i] < min_nnz:
            weights[i] = np.nan
    
    logger.info('Chromosome {0} .. Done'.format(chrom))
    
    return weights, indices

def split_chromosome(clr, chrom, window_size):

    chromsize = clr.chromsizes[chrom]
    res = clr.binsize
    queue = []
    for s in range(0, chromsize, window_size//10*9):
        if chromsize - s > window_size//2*3:
            e = s + window_size
            queue.append((s//res, e//res))
        else:
            if chromsize % res == 0:
                queue.append((s//res, chromsize//res))
            else:
                queue.append((s//res, chromsize//res+1))
            break
    
    return queue

def extract_valid_regions(clr, included_bins, window_size):

    pieces = np.split(included_bins, np.where(np.diff(included_bins)!=1)[0]+1)
    res = clr.binsize
    queue = []
    for arr in pieces:
        start = arr[0] * res
        end = arr[-1] * res
        for s in range(start, end, window_size//10*9):
            if end - s > window_size//2*3:
                e = s + window_size
                queue.append((s//res, e//res))
            else:
                queue.append((s//res, end//res))
                break
    
    return queue

def combine_weights(weights_by_window, ini_weights):

    weights = np.zeros_like(ini_weights)
    for s, e in sorted(weights_by_window):
        for i in range(s, e):
            if weights[i] == 0:
                weights[i] = weights_by_window[(s,e)][0][i-s]
            elif weights[i] is np.nan:
                weights[i] = weights_by_window[(s,e)][0][i-s]
            else:
                weights[i] = (weights[i] + weights_by_window[(s,e)][0][i-s]) / 2
    
    return weights

def write_weights(bias, cool_path, group_path, column_name):

    with h5py.File(cool_path, 'r+') as h5:
        grp = h5[group_path]
        if column_name in grp['bins']:
            del grp['bins'][column_name]
        # add the bias column to the file
        h5opts = dict(compression='gzip', compression_opts=6)
        grp['bins'].create_dataset(column_name, data=bias, **h5opts)

def start_logger(logfil):

     ## Root Logger Configuration
    logger = logging.getLogger()
    logger.setLevel(10)
    console = logging.StreamHandler()
    filehandler = logging.FileHandler(logfil)
    # Set level for Handlers
    console.setLevel('INFO')
    filehandler.setLevel('INFO')
    # Customizing Formatter
    formatter = logging.Formatter(fmt = '%(name)-25s %(levelname)-7s @ %(asctime)s: %(message)s',
                                datefmt = '%m/%d/%y %H:%M:%S')
    
    ## Unified Formatter
    console.setFormatter(formatter)
    filehandler.setFormatter(formatter)
    # Add Handlers
    logger.addHandler(console)
    logger.addHandler(filehandler)

    return logger

def load_BED(infil, res):

    D = defaultdict(set)
    with open(infil, 'r') as source:
        for line in source:
            c, s, e = line.rstrip().split()[:3]
            s, e = int(s), int(e)
            bins = set(range(s//res, (e+res-1)//res))
            D[c].update(bins)
    
    for c in D:
        D[c] = sorted(D[c])
    
    return D