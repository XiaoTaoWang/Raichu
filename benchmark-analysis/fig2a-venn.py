import matplotlib
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.spatial.distance import cdist, euclidean

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def parse_loops(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            c1, s1, e1, c2, s2, e2 = parse[:6]
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            p1 = (s1 + e1) // 2
            p2 = (s2 + e2) // 2
            if not c1 in D:
                D[c1] = []
            D[c1].append((s1, e1, s2, e2, p1, p2))
    
    for c1 in D:
        D[c1] = np.r_[D[c1]]

    return D

def set2dict(S):

    D = {}
    for c, s1, e1, s2, e2 in S:
        if not c in D:
            D[c] = set()
        D[c].add((s1, e1, s2, e2))

    return D

def remove_redundant(obj_common, weight_common):

    obj_common = set2dict(obj_common)
    weight_common = set2dict(weight_common)

    for c in weight_common:
        for w in weight_common[c]:
            w_p1 = (w[0] + w[1]) // 2
            w_p2 = (w[2] + w[3]) // 2
            bychrom = obj_common[c]
            sort_table = []
            for o in bychrom:
                o_p1 = (o[0] + o[1]) // 2
                o_p2 = (o[2] + o[3]) // 2
                d = euclidean([w_p1, w_p2], [o_p1, o_p2])
                sort_table.append((d, o))
            sort_table.sort()
            obj_common[c].discard(sort_table[0][1])
    
    additional_obj = set()
    for c in obj_common:
        for s1, e1, s2, e2 in obj_common[c]:
            additional_obj.add((c, s1, e1, s2, e2))
    
    return additional_obj

def write_loop(loops, outfil):

    with open(outfil, 'w') as out:
        for tmp in loops:
            line = list(map(str, tmp))
            out.write('\t'.join(line)+'\n')

loop_fil1 = 'GM12878_mustache.weight.merge.bedpe'
loop_fil2 = 'GM12878_mustache.obj_weight.merge.bedpe'
dis_cutoff = 50000
ratio_cutoff = 0.2
set1_specific = set()
set2_specific = set()
set1_common = set()
set2_common = set()

loop_set1 = parse_loops(loop_fil1)
loop_set2 = parse_loops(loop_fil2)
# iterate the first loop set
for c in loop_set1:
    for coord in loop_set1[c]:
        if not c in loop_set2:
            set1_specific.add((c,)+tuple(coord[:4]))
        else:
            dis = cdist([coord[[4,5]]], loop_set2[c][:,[4,5]])
            min_dis = dis.min()
            if min_dis < min(ratio_cutoff*abs(coord[5]-coord[4]), dis_cutoff):
                set1_common.add((c,)+tuple(coord[:4]))
            else:
                set1_specific.add((c,)+tuple(coord[:4]))

# iterate the second loop set
for c in loop_set2:
    for coord in loop_set2[c]:
        if not c in loop_set1:
            set2_specific.add((c,)+tuple(coord[:4]))
        else:
            dis = cdist([coord[[4,5]]], loop_set1[c][:,[4,5]])
            min_dis = dis.min()
            if min_dis < min(ratio_cutoff*abs(coord[5]-coord[4]), dis_cutoff):
                set2_common.add((c,)+tuple(coord[:4]))
            else:
                set2_specific.add((c,)+tuple(coord[:4]))


additional = remove_redundant(set2_common, set1_common)
set2_specific = set2_specific | additional

# Venn diagram
loop_set1_total = set1_specific | set1_common
loop_set2_total = set1_common | set2_specific 
fig = plt.figure(figsize=(1.8, 1.8))
ven = venn2([loop_set1_total, loop_set2_total],
            ('ICE', 'Raichu'), set_colors=("#5585C3", "#BBDE87"), alpha=0.4)
for text in ven.set_labels:
    text.set_fontsize(6)
plt.savefig('venn.mustache.ICE-vs-Raichu.svg', dpi=500)

write_loop(set1_specific, 'mustache.ICE-vs-Raichu.ICE-specific.txt')
write_loop(set2_specific, 'mustache.ICE-vs-Raichu.Raichu-specific.txt')
write_loop(set1_common, 'mustache.ICE-vs-Raichu.common.txt')