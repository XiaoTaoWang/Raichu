from HiCLift.liftover import LiftOver

def parse_chip(infil):

    peaks = []
    with open(infil, 'r') as source:
        for i, line in enumerate(source):
            c, s, e = line.rstrip().split()[:3]
            s, e = int(s), int(e)
            peaks.append((c, s, e, i))
    
    return peaks

queue = [
    ('GSE145997_engineered_AA_GATA3_ChIP-seq_pooled.narrowPeak', 'GM12878_AA.GATA3'),
    ('GSE145997_GM12878_CC_GATA3_ChIP-seq_pooled.narrowPeak', 'GM12878_CC.GATA3'),
    ('GSM4349495_engineered_AA_ChIP-seq_H3K27ac.filt.q.0.01.narrowPeak', 'GM12878_AA.H3K27ac'),
    ('GSM4349496_engineered_AA_ChIP-seq_H3K4me1.filt.q.0.01.narrowPeak', 'GM12878_AA.H3K4me1'),
    ('GSM4349497_GM12878_CC_ChIP-seq_H3K27ac.filt.q.0.01.narrowPeak', 'GM12878_CC.H3K27ac'),
    ('GSM4349498_GM12878_CC_ChIP-seq_H3K4me1.filt.q.0.01.narrowPeak', 'GM12878_CC.H3K4me1')
]
lo = LiftOver('hg19', 'hg38')

for infil, outpre in queue:
    peaks = parse_chip(infil)
    outfil = '{0}.peaks.bed'.format(outpre)
    converted = []
    for c, s, e, i in peaks:
        width = e - s
        left = lo.convert_coordinate(c, s)
        right = lo.convert_coordinate(c, e)
        if left is None:
            continue
        if right is None:
            continue
        
        if (len(left) == 1) and (len(right) == 1):
            if (left[0][0] == c) and (right[0][0] == c):
                s_ = left[0][1]
                e_ = right[0][1]
                if e_ - s_ == width:
                    converted.append((c, str(s_), str(e_), str(i)))
    
    with open(outfil, 'w') as out:
        for line in converted:
            out.write('\t'.join(line)+'\n')
