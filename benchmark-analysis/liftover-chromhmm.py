from HiCLift.liftover import LiftOver

def parse_chromhmm(infil):

    data = []
    with open(infil, 'r') as source:
        for line in source:
            if line.startswith('#'):
                continue
            c, s, e, annot = line.rstrip().split()[1:5]
            s, e = int(s), int(e)
            data.append((c, s, e, annot))
    
    return data

data = parse_chromhmm('K562-ChromHMM-hg19.bed')
lo = LiftOver('hg19', 'hg38')
transformed = []
for c, s, e, annot in data:
    s_l = lo.convert_coordinate(c, s)
    e_l = lo.convert_coordinate(c, e)
    if (s_l is None) or (e_l is None):
        continue

    if (len(s_l) != 1) or (len(e_l) != 1):
        continue

    s_l = s_l[0]
    e_l = e_l[0]
    if (s_l[0] != c) or (e_l[0] != c):
        continue

    if (e - s) != (e_l[1] - s_l[1]):
        continue

    tmp = (c, str(s_l[1]), str(e_l[1]), annot)
    transformed.append(tmp)

with open('K562-ChromHMM-hg38.bed', 'w') as out:
    for tmp in transformed:
        out.write('\t'.join(tmp)+'\n')
