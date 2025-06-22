import numpy as np
from HiCLift.liftover import LiftOver
from collections import defaultdict, Counter
from scipy.spatial.distance import cdist
import pickle

def load_loops(infil):

    loops = []
    with open(infil, 'r') as source:
        for line in source:
            c1, s1, e1, s2, e2 = line.strip().split()
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            loops.append((c1, s1, e1, c1, s2, e2))
    
    return loops

def load_ref_loops(fil_list):

    loops = defaultdict(list)
    for infil in fil_list:
        with open(infil, 'r') as source:
            for line in source:
                c1, s1, e1, c2, s2, e2 = line.strip().split()
                s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
                p1 = (s1 + e1) // 2
                p2 = (s2 + e2) // 2
                loops[c1].append((p1, p2))
    
    for c in loops:
        loops[c] = np.array(loops[c])
    
    return loops

def check_in(loop, ref_loops, dis_cutoff=80000, ratio_cutoff=0.3):

    # loop is a tuple of (chrom, p1, p2)
    chrom, p1, p2 = loop
    if chrom not in ref_loops:
        return False
    
    cutoff = min(ratio_cutoff*abs(p2-p1), dis_cutoff)
    ref = ref_loops[chrom]
    dis = cdist([[p1, p2]], ref).ravel().min()
    if dis < cutoff:
        return True
    
    return False

lo = LiftOver('mm10ToHg38.over.chain.gz')
dis_cutoff = 100000
ratio_cutoff = 0.4
n_points = 5001
queue = [
    'ICE-vs-Raichu.ICE-specific.txt',
    'ICE-vs-Raichu.common.txt',
    'ICE-vs-Raichu.Raichu-specific.txt'
]
ref_loops = load_ref_loops(['hNPC.weight.min16.bedpe', 'hNPC.obj_weight_200.min16.bedpe'])
stats = {}
loop_map = {}
for q in queue:
    key = q.split('.')[1]
    loops = load_loops(q)
    n_total = len(loops)
    n_pass = 0
    for c1, s1, e1, c2, s2, e2 in loops:
        pool1 = []
        chroms1 = []
        for p in np.linspace(s1, e1, n_points):
            tmp = lo.convert_coordinate(c1, int(p))
            if (not tmp is None) and (len(tmp) == 1):
                chroms1.append(tmp[0][0])
                pool1.append(tmp[0][1])
        
        if not len(chroms1):
            continue
 
        pool2 = []
        chroms2 = []
        for p in np.linspace(s2, e2, n_points):
            tmp = lo.convert_coordinate(c2, int(p))
            if (not tmp is None) and (len(tmp) == 1):
                pool2.append(tmp[0][1])
                chroms2.append(tmp[0][0])
        
        if not len(chroms2):
            continue

        chrom1 = Counter(chroms1).most_common(1)[0][0]
        chrom2 = Counter(chroms2).most_common(1)[0][0]
        if chrom1 != chrom2:
            continue

        pool1_filtered = [p for p, c in zip(pool1, chroms1) if c == chrom1]
        pool2_filtered = [p for p, c in zip(pool2, chroms2) if c == chrom2]
        p1 = np.mean(pool1_filtered)
        p2 = np.mean(pool2_filtered)
        if p1 > p2:
            p1, p2 = p2, p1
        
        if ((p2 - p1) / (s2 - s1) > 0.8) and ((p2 - p1) / (s2 - s1) < 1.25):
            loop_map[(c1, s1, e1, c2, s2, e2)] = [(chrom1, p1, p1+5000, chrom1, p2, p2+5000), 0]
            if check_in((chrom1, p1, p2), ref_loops, dis_cutoff, ratio_cutoff):
                n_pass += 1
                loop_map[(c1, s1, e1, c2, s2, e2)] = [(chrom1, p1, p1+5000, chrom1, p2, p2+5000), 1]

    stats[key] = (n_total, n_pass)

with open('mNPC-to-hNPC.conservation.stats.loose.pkl', 'wb') as outf:
    pickle.dump([loop_map, stats], outf)
