from collections import Counter
import numpy as np

def parse_bedpe(fil):

    pos1 = {}
    pos2 = {}
    loops = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0]
            if chrom in ['chrM']:
                continue
            c1, s1, e1, c2, s2, e2 = parse[:6]
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            p1 = (s1 + e1) // 2
            p2 = (s2 + e2) // 2
            if not chrom in pos1:
                pos1[chrom] = []
                pos2[chrom] = []
                loops[chrom] = []
            pos1[chrom].append(p1)
            pos2[chrom].append(p2)
            loops[chrom].append((c1, s1, e1, c2, s2, e2))
    
    for c in pos1:
        pos1[c] = np.r_[pos1[c]]
        pos2[c] = np.r_[pos2[c]]
    
    return pos1, pos2, loops

def check_in(p1, p2, ref1, ref2, ref_loops, c, thre=15000):

    hit = []
    if not c in ref1:
        return hit
    tmp1 = ref1[c]
    tmp2 = ref2[c]
    dis1 = np.abs(p1 - tmp1)
    dis2 = np.abs(p2 - tmp2)
    thre = min(thre, 0.2*abs(p2 - p1))
    mask = (dis1 < thre) & (dis2 < thre)
    idx = np.where(mask)[0]
    if idx.size > 0:
        hit = [ref_loops[c][i] for i in idx]
    
    return hit

def read_peakachu(infil):

    D = {}
    with open(infil, 'r') as source:
        for line in source:
            c1, s1, e1, c2, s2, e2, prob = line.rstrip().split()[:7]
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            prob = float(prob)
            D[(c1, s1, s2)] = prob
    
    return D


stringent_list = [
    'GM12878_AA-None-allReps-filtered.obj_weight.hiccups.bedpe',
    'GM12878_CC-None-allReps-filtered.obj_weight.hiccups.bedpe',
    'GM12878_AA-None-allReps-filtered.weight.hiccups.bedpe',
    'GM12878_CC-None-allReps-filtered.weight.hiccups.bedpe'
]

loose_list = [
    'GM12878_AA.obj_weight.hiccups.loose.bedpe',
    'GM12878_CC.obj_weight.hiccups.loose.bedpe',
    'GM12878_AA.weight.hiccups.loose.bedpe',
    'GM12878_CC.weight.hiccups.loose.bedpe'
]

names = ['AA_obj', 'CC_obj', 'AA_weight', 'CC_weight']
outfil = 'union-loops.bedpe'
# load loops
stringent_sets = [parse_bedpe(f) for f in stringent_list]
loose_sets = [parse_bedpe(f) for f in loose_list]
cache = set()
union_loops = [] # 8 columns
for i in range(len(stringent_sets)):
    query1, query2, query_loops = stringent_sets[i]
    for c in query1:
        for p1, p2, tmp_loop in zip(query1[c], query2[c], query_loops[c]):
            if tmp_loop in cache: # current loop have been added in previous steps
                continue
            ID_overlap = [0, 0, 0, 0]
            ID_overlap[i] = 1
            for j in range(len(stringent_sets)):
                if i==j:
                    continue
                ref1, ref2, ref_loops = stringent_sets[j]
                hit = check_in(p1, p2, ref1, ref2, ref_loops, c, thre=50000)
                if len(hit):
                    ID_overlap[j] = 1
                for h_loop in hit:
                    cache.add(h_loop)
                
                # check if the loop was called at more loose parameter settings
                ref1, ref2, ref_loops = loose_sets[j]
                hit = check_in(p1, p2, ref1, ref2, ref_loops, c, thre=50000)
                if len(hit):
                    ID_overlap[j] = 1
            
            overlap_label = ','.join([names[ii] for ii in range(len(ID_overlap)) if ID_overlap[ii]==1])
            union_loops.append(tmp_loop + (overlap_label,))
            cache.add(tmp_loop)

with open(outfil, 'w') as out:
    for loop in sorted(union_loops):
        out.write('\t'.join(list(map(str, loop)))+'\n')
