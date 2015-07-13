import apa
import pybio

def genes_count(species):
    genes = {}
    f = open(apa.path.polyadb_filename(species, filetype="ann"), "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        if data["gid"]!="":
            genes[data["gid"]] = genes.setdefault(data["gid"], 0) + 1
        r = f.readline()
    L = [(c,g) for g,c in genes.items()]
    L.sort(reverse=True)
    return L
    
for species in ["hg19", "mm10", "dm5"]:
    L = genes_count(species)
    counts = {}
    all_genes = 0
    more_10 = 0
    for (c,g) in L:
        all_genes += 1
        if c>10:
            more_10 += 1
        else:
            counts[c] = counts.setdefault(c, 0) + 1
    counts[">10"] = more_10
    print species
    L = [(c, g) for c,g in counts.items()]
    L.sort()
    for count, num_genes in L:
        print "%s\t%s" % (count, num_genes)
    print

