import apa
import os
import sys
import pybio
from orangecontrib.bio import go

def go_enrichment(species, comps_id, pair_type):

    def read_tab_reg(comps_id, pair_type, site_type, reg_type):
        res = []
        reg_type = {"r":"neg", "e":"pos"}[reg_type]
        fname = os.path.join(apa.path.comps_folder, comps_id, "rnamap", "clip0_heat.%s.%s_%s.tab" % (pair_type, site_type, reg_type))
        f = open(fname, "rt")
        header = f.readline().replace("\n", "").replace("\r", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\n", "").replace("\r", "").split("\t")
            data = dict(zip(header, r))
            if int(data["clip"])>0:
                res.append(data["gene_name"])
            r = f.readline()
        f.close()
        return res

    data_file = os.path.join(apa.path.comps_folder, comps_id, "rnamap", "data_%s.tab" % pair_type)

    r_list = []
    e_list = []
    c_list = [] # reference list

    # read in the gene lists
    f = open(data_file, "rt")
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "")
        if r.startswith("genes_repressed="):
            r = r.split("genes_repressed=")[1]
            r_list = eval(r)
        if r.startswith("genes_enhanced="):
            r = r.split("genes_enhanced=")[1]
            e_list = eval(r)
        if r.startswith("genes_controls_down="):
            r = r.split("genes_controls_down=")[1]
            c_list = c_list + eval(r)
        if r.startswith("genes_controls_up="):
            r = r.split("genes_controls_up=")[1]
            c_list = c_list + eval(r)
        r = f.readline()
    f.close()

    c_list = list(set(r_list+e_list+c_list)) # construct reference list from all present genes

    # read in top regulated (bound) genes
    r_proximal = read_tab_reg(comps_id, pair_type, "proximal", "r")
    e_proximal = read_tab_reg(comps_id, pair_type, "proximal", "e")
    r_distal = read_tab_reg(comps_id, pair_type, "distal", "r")
    e_distal = read_tab_reg(comps_id, pair_type, "distal", "e")

    data_file = os.path.join(apa.path.comps_folder, comps_id, "rnamap", "clip0_heat.%s.proximal_neg.png" % pair_type)

    # compute GO enrichment
    ontology = go.Ontology()
    annotations = go.Annotations({"hg19":"hsa", "mm10":"mmu", "mm9":"mmu"}[species], ontology=ontology)

    comps = []
    comps.append((r_list, c_list, "enhanced"))
    comps.append((e_list, c_list, "repressed"))
    comps.append((r_proximal, c_list, "repressed_proximal"))
    comps.append((r_distal, c_list, "repressed_distal"))
    comps.append((e_proximal, c_list, "enhanced_proximal"))
    comps.append((e_distal, c_list, "enhanced_distal"))

    for (gene_set, reference_set, name) in comps:
        print "%s GO analysis: genes=%s, controls=%s" % (name, len(gene_set), len(reference_set))
        res = annotations.get_enriched_terms(gene_set, reference=reference_set, use_fdr=False)
        results = []
        for go_id, (genes, p_value, ref) in res.items():
            go_depth = ontology.term_depth(go_id)
            if go_depth>=4: # ignore lower levels
                results.append((p_value, go_id, go_depth, ontology[go_id].name, ",".join(genes)))

        fname = os.path.join(apa.path.comps_folder, comps_id, "rnamap", "go_%s.tab" % name)
        print "\t-> %s" % fname
        f = open(fname, "wt")
        f.write("#query (%s genes) = " % len(gene_set) + ",".join(gene_set) + "\n")
        f.write("#reference (%s genes) = " % len(reference_set) + ",".join(reference_set) + "\n\n")
        f.write("GO id\tGO term\tGO depth\tp-value\tenriched genes\n")
        results.sort()
        count = 0
        for (p_value, go_id, go_depth, go_term, genes) in results:
            if p_value<0.05:
                count += 1
                f.write("%s\t%s\t%s\t%.4f\t%s\n" % (go_id, go_term, go_depth, p_value, genes))
        f.close()
        print "\t-> %s GO terms enriched" % count
        print
