import apa
import os
import sys
import pybio
import json
#from orangecontrib.bio import go

def is_sorted(l):
    return all(l[i] <= l[i+1] for i in xrange(len(l)-1))

def FDR(p_values, dependent=False, m=None, ordered=False):
    """
    If the user is sure that pvalues as already sorted nondescendingly
    setting ordered=True will make the computation faster.
    """

    if not ordered:
        ordered = is_sorted(p_values)

    if not ordered:
        joined = [ (v,i) for i,v in enumerate(p_values) ]
        joined.sort()
        p_values = [ p[0] for p in joined ]
        indices = [ p[1] for p in joined ]

    if not m:
        m = len(p_values)
    if m <= 0 or not p_values:
        return []

    if dependent: # correct q for dependent tests
        k = c[m-1] if m <= len(c) else math.log(m) + 0.57721566490153286060651209008240243104215933593992
        m = m * k

    tmp_fdrs = [p*m/(i+1.0) for (i, p) in enumerate(p_values)]
    fdrs = []
    cmin = tmp_fdrs[-1]
    for f in reversed(tmp_fdrs):
        cmin = min(f, cmin)
        fdrs.append( cmin)
    fdrs.reverse()

    if not ordered:
        new = [ None ] * len(fdrs)
        for v,i in zip(fdrs, indices):
            new[i] = v
        fdrs = new

    return fdrs

def max_overlap(my_set, sets):
    my_set = set(my_set)
    if len(my_set)==0:
        return 0
    score = 0
    for s in sets:
        i = my_set & s
        score = max(score, (len(i)/float(len(my_set))))
    return int(round(score*100))

def go_enrichment(comps_id):

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

    # get comps and species
    comps = apa.comps.read_comps(comps_id)
    species = comps.species

    for aspect in ["P", "C"]:
        for pair_type in ["tandem", "composite", "skipped"]:

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
                print("{pair_type}:{name} GO analysis: genes={num_genes}, controls={num_control}".format(pair_type=pair_type, name=name, num_genes=len(gene_set), num_control=len(reference_set)))
                res = annotations.get_enriched_terms(gene_set, reference=reference_set, use_fdr=False, aspect=aspect)
                results = []
                for go_id, (genes, p_value, go_genes) in res.items():
                    go_depth = ontology.term_depth(go_id)
                    if go_depth>=4:
                        results.append((p_value, go_id, go_depth, ontology[go_id].name, ",".join(genes), go_genes))
                results.sort()

                fdr_corrected = FDR([rec[0] for rec in results])
                temp_results = []
                for rec, fdr_val in zip(results, fdr_corrected):
                    rec = list(rec)
                    rec.append(fdr_val)
                    temp_results.append(rec)

                # determine regundancy
                considered_sets = []
                results = []
                for (p_value, go_id, go_depth, go_term, genes, go_genes, fdr) in temp_results:
                    redundancy = max_overlap(genes.split(","), considered_sets)
                    results.append((p_value, go_id, go_depth, go_term, genes, go_genes, fdr, redundancy))
                    considered_sets.append(set(genes.split(",")))

                # write JSON file with results
                fname = os.path.join(apa.path.comps_folder, comps_id, "rnamap", "go_%s_%s_%s.json" % (name, pair_type, aspect))
                f = open(fname, "wt")
                f.write(json.dumps(results))
                f.close()

                fname = os.path.join(apa.path.comps_folder, comps_id, "rnamap", "go_%s_%s_%s.tab" % (name, pair_type, aspect))
                print("\t-> %s" % fname)
                f = open(fname, "wt")
                f.write("#query (%s genes) = " % len(gene_set) + ",".join(gene_set) + "\n")
                f.write("#reference (%s genes) = " % len(reference_set) + ",".join(reference_set) + "\n")
                f.write("#aspect = %s\n" % aspect)
                f.write("GO_id\tGO_term\tGO_depth\tGO_genes\tp_value\tfdr\tRedundancy\tenriched_genes\n")
                for (p_value, go_id, go_depth, go_term, genes, go_genes, fdr, redundancy) in results:
                    f.write("%s\t%s\t%s\t%s\t%.4f\t%.4f\t%s\t%s\n" % (go_id, go_term, go_depth, go_genes, p_value, fdr, redundancy, genes))
                f.close()
                print
