import cooler, sys, joblib
import numpy as np
from scipy import sparse
from hicpeaks.callers import *

def hiccups(M, cM, B1, B2, IR, chromLen, Diags, cDiags, num, chrom, pw=[2], ww=[5], maxww=5,
            sig=0.1, sumq=0.01, double_fold=1.75, single_fold=2, maxapart=3000000,
            res=10000, min_marginal_peaks=3, onlyanchor=True, min_local_reads=25):

    # more codes for lower memory
    # use reference instead of creating new arrays
    extDiags_ref = []
    for i in range(num):
        OneDArray = Diags[i]
        extODA = np.zeros(chromLen - i + maxww*2)
        extODA[maxww:-maxww] = OneDArray
        extDiags_ref.append(extODA)
    
    extDiags = {maxww: extDiags_ref}
    for w in range(min(ww), maxww):
        temp = []
        for i in range(num):
            delta = maxww-w
            extODA = extDiags_ref[i][delta:-delta]
            temp.append(extODA)
        extDiags[w] = temp
    
    EDiags = []
    x = np.r_[sorted(IR)]
    for i in x:
        OneDArray = np.ones(chromLen - i) * IR[i]
        EDiags.append(OneDArray)
    
    EM = sparse.diags(EDiags, x, format = 'csr')

    extCDiags_ref = []
    extEDiags_ref = []
    for i in range(x.size):
        extODA_C = np.zeros(chromLen - x[i] + maxww*2)
        extODA_C[maxww:-maxww] = cDiags[i]
        extCDiags_ref.append(extODA_C)
        extODA_E = np.zeros(chromLen - x[i] + maxww*2)
        extODA_E[maxww:-maxww] = EDiags[i]
        extEDiags_ref.append(extODA_E)
    
    extCDiags = {maxww: extCDiags_ref}
    extEDiags = {maxww: extEDiags_ref}
    for w in range(min(ww), maxww):
        tempC = []
        tempE = []
        for i in range(x.size):
            delta = maxww - w
            extODA_C = extCDiags_ref[i][delta:-delta]
            tempC.append(extODA_C)
            extODA_E = extEDiags_ref[i][delta:-delta]
            tempE.append(extODA_E)
        extCDiags[w] = tempC
        extEDiags[w] = tempE
    
    p_w = pw_ww_pairs(pw, ww, maxww)
                
    ## Peak Calling ...    
    vxi, vyi = cM.nonzero()
    Mask = ((vyi - vxi) >= min(ww)) & ((vyi - vxi) <= (maxapart // res))
    vxi = vxi[Mask]
    vyi = vyi[Mask]
    # Here the key indicates the color the original paper uses for corresponding backgrounds
    flocals = ['K', 'Y'] # Order is important
    bSV = {}; bEV = {}
    for pi in pw: # support multiple pw and ww values
        bSV[pi] = {}; bEV[pi] = {}
        for fl in flocals:
            bSV[pi][fl] = np.zeros(vxi.size)
            bEV[pi][fl] = np.zeros(vxi.size)
    
    print('Chrom:{0}, Observed Contact Number: {1}'.format(chrom, vxi.size))

    RefIdx = {}; iniNum = {}
    for pi in pw:
        RefIdx[pi] = np.arange(vxi.size)
        iniNum[pi] = vxi.size
    
    totalNum = vxi.size
    
    print('Chrom:{0}, Two local neighborhoods, two expected matrices ...'.format(chrom))
    bS = {}; bE = {}
    for fl in flocals:
        bS[fl] = sparse.csr_matrix((chromLen, chromLen))
        bE[fl] = sparse.csr_matrix((chromLen, chromLen))
    Reads = sparse.csr_matrix((chromLen, chromLen))
    limitCompute = False
    last_pi = last_wi = 0
    frozen_w = maxww
    for pi, wi in p_w:
        if wi > frozen_w:
            continue
        ps = 2 * pi + 1
        ws = 2 * wi + 1
        print('Chrom:{0},    Peak width:{1}, Donut width:{2}'.format(chrom, pi, wi))
        P1 = set([(i,j) for i in range(wi-pi, ps+wi-pi) for j in range(wi-pi, ps+wi-pi)]) # Center Peak Region
        P_1 = set([(i,j) for i in range(wi+1, ws) for j in range(wi)])
        P_2 = set([(i,j) for i in range(wi+1, ps+wi-pi) for j in range(wi-pi, wi)])
        P2 = P_1 - P_2 # Lower-left Region

        ss = range(ws)
        Pool_Diags = {}
        Pool_EDiags = {}
        Pool_cDiags = {}
        for i in ss:
            for j in ss:
                bgloc = max(abs(i-wi), abs(j-wi)) # mark the radial location on background matrix
                if limitCompute:
                    if ((bgloc<=last_wi) and (bgloc>max(pi,last_pi))) or (bgloc<=min(pi,last_pi)):
                        continue
                Pool_Diags[(i,j)] = []
                Pool_EDiags[(i,j)] = []
                Pool_cDiags[(i,j)] = []

                for oi in range(num):
                    if oi + i - j >= 0:
                        starti = i
                        endi = i + chromLen - (oi + i - j)
                    else:
                        starti = i - (oi + i - j)
                        endi = starti + chromLen + (oi + i - j)
                    Pool_Diags[(i,j)].append(extDiags[wi][oi][starti:endi])
                for oi in range(x.size):
                    if x[oi] + i - j >= 0:
                        starti = i
                        endi = i + chromLen - (x[oi] + i - j)
                    else:
                        starti = i - (x[oi] + i - j)
                        endi = starti + chromLen + (x[oi] + i - j)
                    Pool_EDiags[(i,j)].append(extEDiags[wi][oi][starti:endi])
                    Pool_cDiags[(i,j)].append(extCDiags[wi][oi][starti:endi])

        for key in Pool_Diags:
            bgloc = max(abs(key[0]-wi), abs(key[1]-wi))
            cDiags_matrix = sparse.diags(Pool_cDiags[key], x + (key[0] - key[1]), format = 'csr')
            EDiags_matrix = sparse.diags(Pool_EDiags[key], x + (key[0] - key[1]), format = 'csr')
            if (key[0] != wi) and (key[1] != wi) and (key not in P1) and (key not in P2):
                if (not limitCompute) or (limitCompute and bgloc>last_wi) or (limitCompute and bgloc>pi and bgloc<=last_pi):
                    bS['K'] = bS['K'] + cDiags_matrix
                    bE['K'] = bE['K'] + EDiags_matrix
                else:
                    bS['K'] = bS['K'] - cDiags_matrix
                    bE['K'] = bE['K'] - EDiags_matrix
            if key in P2:
                if (not limitCompute) or (limitCompute and bgloc>last_wi) or (limitCompute and bgloc>pi and bgloc<=last_pi):
                    bS['K'] = bS['K'] + cDiags_matrix
                    bE['K'] = bE['K'] + EDiags_matrix
                    bS['Y'] = bS['Y'] + cDiags_matrix
                    bE['Y'] = bE['Y'] + EDiags_matrix
                else:
                    bS['K'] = bS['K'] - cDiags_matrix
                    bE['K'] = bE['K'] - EDiags_matrix
                    bS['Y'] = bS['Y'] - cDiags_matrix
                    bE['Y'] = bE['Y'] - EDiags_matrix
                if (not limitCompute) or (limitCompute and pi==min(pw) and bgloc>last_wi):
                    Reads = Reads + sparse.diags(Pool_Diags[key], np.arange(num) + (key[0] - key[1]), format = 'csr')
            
        limitCompute = True
        last_pi, last_wi = pi, wi
                
        Txi = vxi[RefIdx[pi]]
        Tyi = vyi[RefIdx[pi]]
        RNums = np.array(Reads[Txi, Tyi]).ravel()
        EIdx = RefIdx[pi][RNums >= min_local_reads]
        print('Chrom:{0},    ({1},{2}) Valid Contact Number from This Loop: {3}'.format(chrom, pi, wi, EIdx.size))
        Valid_Ratio = EIdx.size/float(iniNum[pi])
        Exi = vxi[EIdx]
        Eyi = vyi[EIdx]
        for fl in flocals:
            bSV[pi][fl][EIdx] = np.array(bS[fl][Exi, Eyi]).ravel()
            bEV[pi][fl][EIdx] = np.array(bE[fl][Exi, Eyi]).ravel()
                
        RefIdx[pi] = RefIdx[pi][RNums < min_local_reads]
            
        iniNum[pi] = RefIdx[pi].size

        left_Ratio = iniNum[pi]/float(totalNum)

        print('Chrom:{0},    ({1},{2}) Total Valid Ratio after This Loop: {3:.3f}'.format(chrom, pi, wi, 1-left_Ratio))
        
        if (Valid_Ratio < 0.3) and (wi >= max(ww)):
            print('Chrom:{0},    ({1},{2}) Ratio of valid contact is too small, assign maximum donut width to {3} ...'.format(chrom, pi, wi, wi))
            frozen_w = wi
        
        if (left_Ratio < 0.03) and (wi >= max(ww)):
            print('Chrom:{0},    ({1},{2}) Very few or no contacts are left, assign maximum donut width to {3} ...'.format(chrom, pi, wi, wi))
            frozen_w = wi
        
        if wi<frozen_w:
            print('Chrom:{0},    ({1},{2}) {3} Contacts will get into next loop ...'.format(chrom, pi, wi, RefIdx[pi].size))
    
    pixel_table = {} # Store combined peak list
    
    print('Chrom:{0}, Poisson Models and Benjamini-Hochberg Correcting for lambda chunks ...'.format(chrom))
    Description = {'K': 'Donut backgrounds', 'Y': 'Lower-left backgrounds'}
    gaps = set(np.where(np.array(cM.sum(axis=1)).ravel() == 0)[0])
    for pi, wi in zip(pw, ww):
        xpos = {}; ypos = {}; Ovalues = {}; ICE = {}
        Fold = {}; pvalues = {}; qvalues = {}
        for fl in flocals:
            print('Chrom:{0},    Peak width:{1}, Donut width:{2}, {3} ...'.format(chrom, pi, wi, Description[fl]))
            Mask = (bEV[pi][fl] != 0) & (vyi - vxi >= wi)
            tmp = sparse.lil_matrix((chromLen, chromLen))
            tmp[vxi[Mask],vyi[Mask]] = bSV[pi][fl][Mask] / bEV[pi][fl][Mask]
            cEM = EM.multiply(tmp.tocsr())
            xi, yi = cEM.nonzero()
            Evalues = np.array(cEM[xi, yi]).ravel() * B1[xi] * B2[yi]
            Mask = Evalues > 0
            Evalues = Evalues[Mask]
            xi = xi[Mask]
            yi = yi[Mask]
            Ovalues[fl] = np.array(M[xi, yi]).ravel()
            ICE[fl] = np.array(cM[xi, yi]).ravel()
            Fold[fl] =  Ovalues[fl] / Evalues
            print('Chrom:{0},    ({1},{2}), Valid contact number: {3}'.format(chrom, pi, wi, xi.size))
        
            pvalue = np.ones(xi.size)
            qvalue = np.ones(xi.size)
        
            print('Chrom:{0},    ({1},{2}), Lambda chunking ...'.format(chrom, pi, wi))
            Bvalues = B1[xi] * B2[yi]
            chunks = lambdachunk(Bvalues)
            print('Chrom:{0},    ({1},{2}), Number of chunks: {3}'.format(chrom, pi, wi, len(chunks)))
            for chunk in chunks:
                print('Chrom:{0},        ({1},{2}), lv: {3:.3g}, rv: {4:.3g}, Num: {5}'.format(chrom, pi, wi, chunk[0], chunk[1], chunk[2].size))
                if chunk[2].size > 0:
                    Poiss = poisson(Evalues[chunk[2]])
                    print('Chrom:{0},        ({1},{2}), Assign P values ...'.format(chrom, pi, wi))
                    chunkP = 1 - Poiss.cdf(Ovalues[fl][chunk[2]])
                    pvalue[chunk[2]] = chunkP
                    print('Chrom:{0},        ({1},{2}), Multiple testing ...'.format(chrom, pi, wi))
                    cResults = multipletests(chunkP, alpha = sig, method = 'fdr_bh')
                    cP = cResults[1] # Corrected Pvalue
                    qvalue[chunk[2]] = cP
                else:
                    print('Chrom:{0},        ({1},{2}), Skipping ...'.format(chrom, pi, wi))
        
            reject = qvalue <= sig
            qvalue = qvalue[reject]
            pvalue = pvalue[reject]
            Ovalues[fl] = Ovalues[fl][reject]
            ICE[fl] = ICE[fl][reject]
            Evalues = Evalues[reject]
            Fold[fl] = Fold[fl][reject]
            xi = xi[reject]
            yi = yi[reject]
        
            print('Chrom:{0},    ({1},{2}), Remove Gap Effects ...'.format(chrom, pi, wi))
        
            if len(gaps) > 0:
                fIdx = []
                for i in np.arange(xi.size):
                    lower = (xi[i] - min(ww)) if (xi[i] > min(ww)) else 0
                    upper = (xi[i] + min(ww)) if ((xi[i] + min(ww)) < chromLen) else (chromLen - 1)
                    cregion_1 = range(lower, upper)
                    lower = (yi[i] - min(ww)) if (yi[i] > min(ww)) else 0
                    upper = (yi[i] + min(ww)) if ((yi[i] + min(ww)) < chromLen) else (chromLen - 1)
                    cregion_2 = range(lower, upper)
                    cregion = set(cregion_1) | set(cregion_2)
                    intersect = cregion & gaps
                    if len(intersect) == 0:
                        fIdx.append(i)
        
                xi = xi[fIdx]
                yi = yi[fIdx]
                Ovalues[fl] = Ovalues[fl][fIdx]
                ICE[fl] = ICE[fl][fIdx]
                pvalue = pvalue[fIdx]
                qvalue = qvalue[fIdx]
                Fold[fl] = Fold[fl][fIdx]
                Evalues = Evalues[fIdx]
        
            xpos[fl] = xi
            ypos[fl] = yi
            pvalues[fl] = pvalue
            qvalues[fl] = qvalue
    
        print('Chrom:{0},    Peak width:{1}, Donut width:{2}, Combine two local filters ...'.format(chrom, pi, wi))

        preDonuts = dict(zip(zip(xpos['K'], ypos['K']), zip(Ovalues['K'], ICE['K'], Fold['K'], pvalues['K'], qvalues['K'])))
        preLL = dict(zip(zip(xpos['Y'], ypos['Y']), zip(Ovalues['Y'], ICE['Y'], Fold['Y'], pvalues['Y'], qvalues['Y'])))
    
        commonPos = set(preDonuts.keys()) & set(preLL.keys())
        postcheck = set(preDonuts.keys()) - set(preLL.keys()) # handle special cases for new peak calling
        for ci, cj in postcheck:
            if cEM[ci,cj]==0: # corresponds to lower-left
                commonPos.add((ci,cj))
        
        Donuts = {}; LL = {}
        for ci, cj in commonPos:
            Donuts[(ci,cj)] = preDonuts[(ci,cj)]
            if (ci,cj) in preLL:
                LL[(ci,cj)] = preLL[(ci,cj)]
            else:
                LL[(ci,cj)] = preDonuts[(ci,cj)]

        for pixel in Donuts:
            donut, ll = Donuts[pixel], LL[pixel]
            key = (pixel[0]*res, pixel[1]*res)
            if (donut[1]>double_fold) and (ll[1]>double_fold) and ((donut[1]>single_fold) or (ll[1]>single_fold)):
                if not key in pixel_table:
                    pixel_table[key] = key + (0,) + donut + ll[2:]
                else:
                    if (donut[-1]<pixel_table[key][7]) and (ll[-1]<pixel_table[key][10]):
                        pixel_table[key] = key + (0,) + donut + ll[2:]

    return pixel_table

# parameter lines
uri = 'K562_H3K27ac-MboI-allReps-filtered.mcool::resolutions/5000'
balance = 'obj_weight3'
chrom = 'chr21'
outfil = 'K562-H3K27ac.{0}.{1}.hiccups.txt'.format(chrom, balance)
pw = [2,4]
ww = [5,7]
maxww = 7
maxapart = 3000000

# prepare data for HiCCUPS
clr = cooler.Cooler(uri)
res = clr.binsize
H = clr.matrix(balance=False, sparse=True).fetch(chrom)
cHeatMap = clr.matrix(balance=balance, sparse=True).fetch(chrom)
chromLen = H.shape[0]
num = maxapart // res + maxww + 1
Diags = [H.diagonal(i) for i in np.arange(num)]
M = sparse.diags(Diags, np.arange(num), format='csr')
x = np.arange(min(ww), num)
IR = {}
cDiags = []
for i in x:
    diag = cHeatMap.diagonal(i)
    mask = np.isnan(diag)
    notnan = diag[np.logical_not(mask)]
    IR[i] = notnan.mean()
    diag[mask] = 0
    cDiags.append(diag)
cM = sparse.diags(cDiags, x, format='csr')

joblib.dump(IR, outfil.replace('.txt', '.pkl'))
del H, cHeatMap

tmp = clr.bins().fetch(chrom)[balance].values
mask = np.logical_not((tmp==0) | np.isnan(tmp))
biases = np.zeros_like(tmp)
biases[mask] = 1/tmp[mask]

# run HiCCUPS
pixel_table = hiccups(M, cM, biases, biases, IR, chromLen, Diags, cDiags, num, chrom, pw=pw, ww=ww,
                      maxww=maxww, sig=1, sumq=1, maxapart=maxapart, double_fold=0, single_fold=0,
                      res=res, min_marginal_peaks=3, onlyanchor=True)

with open(outfil, 'w') as out:
    for pixel in pixel_table:
        lineFormat = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.3g}\t{7:.3g}\t{8:.3g}\t{9:.3g}\t{10:.3g}\t{11:.3g}\t{12:.3g}\t{13:.3g}\n'
        tmp = pixel_table[pixel]
        content = (chrom, pixel[0], pixel[0]+res, chrom, pixel[1], pixel[1]+res) + tmp[3:]
        line = lineFormat.format(*content)
        out.write(line)

