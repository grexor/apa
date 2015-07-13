import apa
import pybio

dist = {}
data = pybio.data.Bedgraph(apa.path.polyadb_filename("hg19.tian", filetype="bed"))
for chr, strand_data in data.raw.items():
    for strand, pos_data in strand_data.items():
        L = [pos for pos, cDNA in pos_data.items()]
        L.sort()
        for pos1, pos2 in zip(L, L[1:]):
            distance = pos2-pos1+1
            if distance<=100:
                dist[100] = dist.setdefault(100, 0) + 1
L = [(c,d) for d,c in dist.items()]
L.sort()
for (c,d) in L:
    print "%s of sites <  %snt apart" % (c,d)
