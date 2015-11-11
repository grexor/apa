import apa
import pybio

def update_vectors(av,tv,cv,gv,seq):
    a = []
    c = []
    t = []
    g = []
    for i in range(0, len(seq)):
        if seq[i]=="A":
            a.append(1)
            c.append(0)
            t.append(0)
            g.append(0)
        elif seq[i]=="C":
            a.append(0)
            c.append(1)
            t.append(0)
            g.append(0)
        elif seq[i]=="T":
            a.append(0)
            c.append(0)
            t.append(1)
            g.append(0)
        elif seq[i]=="G":
            a.append(0)
            c.append(0)
            t.append(0)
            g.append(1)
        else:
            a.append(0)
            c.append(0)
            t.append(0)
            g.append(0)
    if av==None:
        av = a
    else:
        av = [x+y for x,y in zip(av, a)]
    if tv==None:
        tv = t
    else:
        tv = [x+y for x,y in zip(tv, t)]
    if cv==None:
        cv = c
    else:
        cv = [x+y for x,y in zip(cv, c)]
    if gv==None:
        gv = g
    else:
        gv = [x+y for x,y in zip(gv, g)]
    return av,tv,cv,gv

def norm_vectors(a,t,c,g, sites):
    a = [x/float(sites) for x in a]
    t = [x/float(sites) for x in t]
    c = [x/float(sites) for x in c]
    g = [x/float(sites) for x in g]
    return a,t,c,g

for species in ["hg19.tian", "hg19", "mm10", "dm5"]:
    sites = 0
    a = None
    c = None
    t = None
    g = None
    data = pybio.data.Bedgraph(apa.path.polyadb_filename(species, filetype="bed"))
    for chr, strand_data in data.raw.items():
        for strand, pos_data in strand_data.items():
            for pos, cDNA in pos_data.items():
                sites += 1
                seq = pybio.genomes.seq(species.rstrip(".tian"), chr, strand, pos, start=-50, stop=50)
                assert(len(seq)==101)
                a,t,c,g = update_vectors(a,t,c,g,seq)
    print "%s, sites=%s" % (species, sites)
    a,t,c,g = norm_vectors(a,t,c,g, sites)
    f = open("%s.atcg.tab" % species, "wt")
    for i in range(0, len(a)):
        f.write("%s\t%s\t%s\t%s\n" % (a[i], t[i], c[i], g[i]))
    f.close()
