import Orange.bio.go as obiGO, obiProb
import misc, sys

file_reference = sys.argv[1]
file_group = sys.argv[2]
our_label = sys.argv[3]

def readExonFormat(fname):
    tmpl = misc.readList(fname)
    retl = []
    g2e_map = {}
    for e in tmpl:
        l = e.split('\t')
        assert(len(l) == 2)
        [ex, g] = l
        tl = g2e_map.get(g, set())
        nid = "%s_%s" % (g, ex)
        tl.add( nid)
        g2e_map[g] = tl
        retl.append(nid)
    return retl, g2e_map

ref1, dic1 = readExonFormat(file_reference)
cat1, _ = readExonFormat(file_group)

print "gene_exon records in %s (%s): %s" % (our_label, len(cat1), ", ".join(cat1))
print "control exp:", len(ref1)

ontology = obiGO.Ontology.Load()
annot1 = obiGO.Annotations.Load("goa_human", ontology=ontology)
annot1.RemapGenes(dic1)
print "loaded"
sys.stdout.flush()

def maxInter(tmpl, tmpll):
    maxI = 0.0
    tmps = set(tmpl)
    if len(tmps) == 0:
        return 0.0
    for tmps2 in tmpll:
        iset = tmps & tmps2
        maxI = max(maxI, float(len(iset))/float(len(tmps)))
    return maxI

summaryByGO = {}

def calcAndReportFor(genes, annotations, catlabel, reference, aspect, minS, maxS):
    res = annotations.get_enriched_terms(genes, reference=reference, aspect=aspect, use_fdr=False) 
    #results = [(p_value, GOId, ref, genes) for (GOId, (genes, p_value, ref)) in res.items()]
    ignore_lower_levels = 4
    results = [(p_value, GOId, ref, genes) for (GOId, (genes, p_value, ref)) in res.items() if ontology.term_depth(GOId) >= ignore_lower_levels]
    results.sort()
    print "number of tests:", len(results)
    pvals = [p_value for (p_value, GOId, ref, genes) in results]
    fdrs = obiProb.FDR(pvals)
    print "p\tfdr\tterm\tredundant\tterm level\tGO ID\treference genes\tmatching cluster genes\tmatching cluster genes"
    glisted = []
    for (fdr, (p_value, GOId, ref, genes)) in zip(fdrs, results):
        genes.sort()
#		if p_value < 0.1 or fdr < 0.01:
        if 1 or (fdr < 0.05 and len(genes) >= 3 and ref >= minS and ref <= maxS):
            mi =  maxInter(genes, glisted)
            if mi > 0.7:
                red = "RED (%2.0f)" % (mi*100.0)
            else:
                red = "%2.0f" % (mi*100.0)
            print "%.3g\t%.3g\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (p_value, fdr, ontology.terms[GOId].name, red, ontology.term_depth(GOId), GOId, ref, len(genes), ", ".join(genes))
            _, _, tmpd = summaryByGO.get(GOId, (aspect, ontology.terms[GOId].name, {}))
            tmpd[catlabel] = (fdr, genes)
            summaryByGO[GOId] = (aspect, ontology.terms[GOId].name, tmpd)
            glisted.append( set(genes))
    print

for (aspect, minS, maxS) in [("P", 50, 700), ("C", 500, 3000)]: #, "F"]:
    print "---------------------------------------------------------"
    print "aspect:", aspect

    print "%s vs reference:" % (our_label)
    calcAndReportFor(cat1, annot1, 'splice_cat1', ref1, aspect, minS, maxS)
    summaryByGO['splice_cat1'] = cat1

    print

misc.writeDictionary("sum_all_splicing.txt", summaryByGO)
