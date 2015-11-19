import matplotlib
matplotlib.use("Agg", warn=False)
import matplotlib.pyplot as plt
import math
from matplotlib import cm as CM
import numpy
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import pybio
import apa
import os
import fisher as fisher_module
import shutil
import json
import glob
import fisher
import sys
import operator
import itertools
from collections import Counter

class Comps:
    def __init__(self, comps_id=None, iCLIP_filename=None):
        self.comps_id = comps_id
        self.test = []
        self.control = []
        self.test_name = ""
        self.control_name = ""
        self.species = ""
        self.iCLIP_filename = iCLIP_filename
        self.cDNA_thr = 5 # at least cDNA
        self.presence_thr = 2.0 # for at least half of experiments
        self.pc_thr = 0.1
        self.fisher_thr = 0.1
        self.control_thr = 0.025
        self.pair_dist = 450
        self.exp_data = {}
        self.polya_db = None
        self.polya_db_filter = ["strong"] # default
        self.deepbind = None
        self.rnamaps = []
        self.db_type="cs" # cs = cleavage site, pas = polyadenylation signal, cs = default

    def __str__(self):
        print "comps_id = %s" % (self.comps_id)
        for test_id, test_data in self.test.items():
            print "%s %s" % (test_id, test_data)
        for control_id, control_data in self.control.items():
            print "%s %s" % (control_id, control_data)
        print "iCLIP = %s" % self.iCLIP_filename
        return ""

def read_comps(comps_id):
    comps = Comps(comps_id)
    config_file = apa.path.comps_config_filename(comps_id)
    if not os.path.exists(config_file):
        return None
    f = open(config_file, "rt")
    r = f.readline()
    while r.startswith("#"):
        r = f.readline()
    header = r.replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "")
        r = r.split(" #")[0] # remove comments
        r = r.rstrip() # remove whitespace characters from end of string
        r = r.split("\t")
        data = dict(zip(header, r))
        if r==[""]:
            r = f.readline()
            continue
        if r[0].startswith("#"):
            r = f.readline()
            continue
        if r[0].startswith("pc_thr:"):
            comps.pc_thr = float(r[0].split("pc_thr:")[1])
            r = f.readline()
            continue
        if r[0].startswith("fisher_thr"):
            comps.fisher_thr = float(r[0].split("fisher_thr:")[1])
            r = f.readline()
            continue
        if r[0].startswith("control_thr"):
            comps.control_thr = float(r[0].split("control_thr:")[1])
            r = f.readline()
            continue
        if r[0].startswith("pair_dist"):
            comps.pair_dist = float(r[0].split("pair_dist:")[1])
            r = f.readline()
            continue
        if r[0].startswith("iCLIP"):
            comps.iCLIP_filename = r[0].split("iCLIP:")[1]
            r = f.readline()
            continue
        if r[0].startswith("cDNA_thr"):
            comps.cDNA_thr = int(r[0].split("cDNA_thr:")[1])
            r = f.readline()
            continue
        if r[0].startswith("presence_thr:"):
            comps.presence_thr = int(r[0].split("presence_thr:")[1])
            r = f.readline()
            continue
        if r[0].startswith("control_name:"):
            comps.control_name = r[0].split("control_name:")[1]
            r = f.readline()
            continue
        if r[0].startswith("test_name:"):
            comps.test_name = r[0].split("test_name:")[1]
            r = f.readline()
            continue
        if r[0].startswith("polya_db:"):
            comps.polya_db = r[0].split("polya_db:")[1]
            r = f.readline()
            continue
        if r[0].startswith("polya_db_filter:"):
            comps.polya_db_filter = eval(r[0].split("polya_db_filter:")[1])
            r = f.readline()
            continue
        if r[0].startswith("deepbind:"):
            comps.deepbind = r[0].split("deepbind:")[1]
            r = f.readline()
            continue
        if r[0].startswith("rnamaps:"):
            comps.rnamaps = r[0].split("rnamaps:")[1].split(",")
            r = f.readline()
            continue
        if r[0].startswith("db_type:"):
            comps.db_type = r[0].split("db_type:")[1]
            r = f.readline()
            continue
        id = data["id"]
        experiments = data["experiments"].split(",")
        name = data["name"]
        if id.startswith("c"):
            comps.control.append((id, experiments, name))
        if id.startswith("t"):
            comps.test.append((id, experiments, name))
        r = f.readline()
    f.close()

    # determine comps species
    species = set()
    all_exp = set()
    for (_, exp, _) in comps.test:
        all_exp.update(exp)
    for (_, exp, _) in comps.control:
        all_exp.update(exp)
    for exp in all_exp:
        lib_id = exp[:exp.rfind("_")]
        exp_id = int(exp.split("_")[-1][1:])
        species.add(apa.annotation.libs[lib_id].experiments[exp_id]["map_to"])
    if len(comps.control)>0 and len(comps.test)>0:
        assert(len(species)==1)
        comps.species = species.pop()

    # finally, load experiment annotation into comps object
    for exp in all_exp:
        lib_id = exp[:exp.rfind("_")]
        exp_id = int(exp.split("_")[-1][1:])
        comps.exp_data[exp] = apa.annotation.libs[lib_id].experiments[exp_id]
    return comps

def save(comps):
    config_file = apa.path.comps_config_filename(comps.comps_id)
    f = open(config_file, "wt")
    f.write("\t".join(["id", "experiments", "name"]) + "\n")
    for (id, exp_list, name) in comps.control:
        f.write("\t".join(str(x) for x in [id, ",".join(exp_list), name]) + "\n")
    for (id, exp_list, name) in comps.test:
        f.write("\t".join(str(x) for x in [id, ",".join(exp_list), name]) + "\n")
    if comps.iCLIP_filename!="":
        f.write("iCLIP:%s" % str(comps.iCLIP_filename))
    f.close()
    return True

def process_comps(comps_id):

    # clean
    assert(len(apa.path.comps_folder)>0)
    files = glob.glob(os.path.join(apa.path.comps_folder, comps_id, "*"))
    for f in files:
        if not f.endswith(".config"):
            print "removing: %s" % f
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)

    comps = read_comps(comps_id) # study data
    pybio.genomes.load(comps.species)

    # if there is a polya-db specified in the comparison, load the positions into the filter
    # (strong, weak, less)
    poly_filter = {}
    if comps.polya_db!=None and comps.db_type=="cs":
        polydb = apa.polya.read(comps.polya_db)

    replicates = []
    expression = {} # keys = c1, c2, c3, t1, t2, t3...items = bedgraph files

    for (comp_id, experiments, comp_name) in comps.control+comps.test:
        key = "%s:%s" % (comp_id, comp_name)
        print "reading: ", key
        expression[comp_id] = pybio.data.Bedgraph()
        replicates.append((comp_id, key)) # (short_name, long_name)
        for id in experiments:
            lib_id = id[:id.rfind("_")]
            exp_id = int(id.split("_")[-1][1:])
            if comps.db_type=="cs":
                e_filename = apa.path.e_filename(lib_id, exp_id)
            elif comps.db_type=="pas":
                e_filename = apa.path.e_filename(lib_id, exp_id, filetype="pas")
            expression[comp_id].load(e_filename)
            print "adding: %s to %s" % (id, comp_id)
        print

    replicates.sort()

    # save bedgraph files for c1, t1, ...
    beds_folder = os.path.join(apa.path.comps_folder, comps_id, "beds")
    if not os.path.exists(beds_folder):
        os.makedirs(beds_folder)
    for (rshort, rlong) in replicates:
        fname = os.path.join(beds_folder, "%s.%s.bed" % (comps_id, rshort))
        expression[rshort].save(fname, track_id="%s.%s" % (comps_id, rshort), genome=comps.species)

    # combine all into single file
    os.system("cat %s/*.bed > %s" % (beds_folder, os.path.join(beds_folder, "%s_all.bed" % comps_id)))

    # finally, save two groups of summed up experimental groups: control (1) and test (2)
    for (sample_name, sample_data) in [("control", comps.control), ("test", comps.test)]:
        b = pybio.data.Bedgraph()
        for (comp_id, experiments, comp_name) in sample_data:
            for id in experiments:
                lib_id = id[:id.rfind("_")]
                exp_id = int(id.split("_")[-1][1:])
                if comps.db_type=="cs":
                    e_filename = apa.path.e_filename(lib_id, exp_id)
                elif comps.db_type=="pas":
                    e_filename = apa.path.e_filename(lib_id, exp_id, filetype="pas")
                b.load(e_filename)
        bed_filename = os.path.join(beds_folder, "%s.%s_all.bed" % (comps_id, sample_name))
        b.save(bed_filename, track_id="%s.%s_all" % (comps_id, sample_name), genome=comps.species)

    # find common list of positions
    # filter positions
    positions = {}
    for id, bg in expression.items():
        for chr, strand_data in bg.raw.items():
            positions.setdefault(chr, {})
            for strand, pos_set in strand_data.items():
                positions.setdefault(chr, {}).setdefault(strand, set())
                valid_positions = set()
                for pos in pos_set:
                    if comps.polya_db!=None and comps.db_type=="cs":
                        polya_db_item = polydb.get((chr, strand, pos), None)
                        if polya_db_item.get("pas_type", None) in comps.polya_db_filter:
                            valid_positions.add(pos)
                    else:
                        valid_positions.add(pos)
                positions[chr][strand] = positions[chr][strand].union(valid_positions)

    # organize polya sites / pas signals inside genes
    gsites = {}
    for chr, strand_data in positions.items():
        for strand, pos_set in strand_data.items():
            pos_set = list(pos_set)
            for pos in pos_set:
                gid_up, gid, gid_down, gid_interval = apa.polya.annotate_position(comps.species, chr, strand, pos)
                if gid==None:
                    continue # only consider polya sites inside genes
                sites = gsites.get(gid, {})
                cDNA_sum = 0
                expression_vector = []
                site_data = {"pos":pos} # store position also in ES
                for (rshort, _) in replicates:
                    bg = expression[rshort]
                    cDNA = bg.get_value(chr, strand, pos)
                    cDNA_sum += cDNA
                    expression_vector.append(cDNA)
                    site_data[rshort] = cDNA

                # filter lowly expressed sites
                expression_vector = [1 if cDNA>=comps.cDNA_thr else 0 for cDNA in expression_vector]
                if sum(expression_vector) < len(expression_vector)/comps.presence_thr:
                    continue

                site_data["cDNA_sum"] = int(cDNA_sum)
                sites[pos] = site_data
                gsites[gid] = sites

    # expression gene level
    fname = apa.path.comps_expression_filename(comps_id)
    f_genes = open(fname, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype"]
    header += ["sites_position_cDNA", "sites_num", "cDNA_sum"]
    for (rshort, rlong) in replicates:
        header.append(rlong)
    f_genes.write("\t".join(header)+"\n")
    for gene_id, sites in gsites.items():
        gene = apa.polya.get_gene(comps.species, gene_id)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)
        row = [chr, strand, gene_locus, gene_id, gene["gene_name"], gene["gene_biotype"]]
        row.append(sites)
        row.append(len(sites))
        row_comp = []
        cDNA_gtotal = 0
        for (rshort, _) in replicates:
            cDNA_rtotal = 0
            for site_pos, es in sites.items():
                val = es.get(rshort, 0)
                cDNA_rtotal += val
            cDNA_gtotal += cDNA_rtotal
            row_comp.append(cDNA_rtotal)
        row.append(cDNA_gtotal)
        row += row_comp
        f_genes.write("\t".join([str(x) for x in row]) + "\n")
    f_genes.close()

    # expression site level
    gene_sites = {}
    fname = apa.path.comps_expression_filename(comps_id, filetype="sites")
    f_sites = open(fname, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "site_pos", "cDNA_sum"]
    for (rshort, rlong) in replicates:
        header.append(rlong)
    f_sites.write("\t".join(header)+"\n")
    for gene_id, sites in gsites.items():
        gene = apa.polya.get_gene(comps.species, gene_id)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)
        for site_pos, es in sites.items():
            row_1 = [chr, strand, gene_locus, gene_id, gene["gene_name"], gene["gene_biotype"], site_pos]
            row_2 = []
            cDNA_sum = 0
            for (rshort, rlong) in replicates:
                val = es.get(rshort, 0)
                row_2.append(val)
                cDNA_sum += val
            row = row_1 + [cDNA_sum] + row_2
            f_sites.write("\t".join([str(x) for x in row]) + "\n")
    f_sites.close()

    # voom analysis
    expression_sites_fname = apa.path.comps_expression_filename(comps_id, filetype="sites")
    R_file = os.path.join(apa.path.root_folder, "comps", "comps_voom.R")
    command = "R --vanilla --args %s %s %s %s %s < %s" % (os.path.join(apa.path.comps_folder, comps_id), comps_id, expression_sites_fname, len(comps.control), len(comps.test), R_file)
    print command
    pybio.utils.Cmd(command).run()

    # differential gene expression with edgeR
    R_file = os.path.join(apa.path.root_folder, "comps", "comps_edgeR.R")
    input_fname = apa.path.comps_expression_filename(comps_id)
    output_fname = os.path.join(apa.path.comps_folder, comps_id, "%s.genes_de.tab" % comps_id)
    command = "R --vanilla --args %s %s %s %s < %s" % (input_fname, output_fname, len(comps.control), len(comps.test), R_file)
    print command
    pybio.utils.Cmd(command).run()

    # replicates cluster analysis
    R_file = os.path.join(apa.path.root_folder, "comps", "comps_cluster.R")
    input_fname = apa.path.comps_expression_filename(comps_id)
    output_fname = os.path.join(apa.path.comps_folder, comps_id, "%s.cluster_genes.svg" % comps_id)
    command = "R --vanilla --args %s %s %s %s '%s (gene level)' < %s" % (input_fname, output_fname, len(comps.control), len(comps.test), comps_id, R_file)
    print command
    pybio.utils.Cmd(command).run()

    # heatmap of differentially expressed genes
    R_file = os.path.join(apa.path.root_folder, "comps", "comps_heatmap_genes.R")
    input_fname = os.path.join(apa.path.comps_folder, comps_id, "%s.genes_de.tab" % comps_id)
    output_fname = os.path.join(apa.path.comps_folder, comps_id, "%s.heatmap_genes" % comps_id)
    command = "R --vanilla --args %s %s %s %s %s < %s" % (input_fname, output_fname, len(comps.control), len(comps.test), comps_id, R_file)
    print command
    pybio.utils.Cmd(command).run()

    # pairs_de file
    pairs_filename = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
    f_pairs = open(pairs_filename, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "num_sites", "siteup_pos", "siteup_exp", "siteup_UG", "sitedown_pos", "sitedown_exp", "sitedown_UG"]
    if comps.control_name!="":
        header.append("up_control [%s]" % comps.control_name)
        header.append("up_control_sum [%s]" % comps.control_name)
        header.append("down_control [%s]" % comps.control_name)
        header.append("down_control_sum [%s]" % comps.control_name)
    else:
        header.append("up_control")
        header.append("up_control_sum")
        header.append("down_control")
        header.append("down_control_sum")

    if comps.test_name!="":
        header.append("up_test [%s]" % comps.test_name)
        header.append("up_test_sum [%s]" % comps.test_name)
        header.append("down_test [%s]" % comps.test_name)
        header.append("down_test_sum [%s]" % comps.test_name)
    else:
        header.append("up.test")
        header.append("up_test_sum")
        header.append("down.test")
        header.append("down_test_sum")

    header += ["pc", "fisher", "pair_type"]

    f_pairs.write("\t".join(header)+"\n")
    results = []
    for gene_id, sites in gsites.items():
        gene = apa.polya.get_gene(comps.species, gene_id)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)

        considered_types = set()

        if len(sites)>=2:
            L = [(sites[pos]["cDNA_sum"], sites[pos]) for pos in sites.keys()]
            L.sort(reverse=True) # sort by control+test expression
            # todo1: here implement site selection
            site_pairs = []
            site_pairs_stats = Counter()
            for (_, site1), (_, site2) in list(itertools.combinations(L, 2)): # skip the cDNA sum used to sort the list
                #print site1["pos"], site2["pos"]
                pair_type = apa.polya.annotate_pair(comps.species, chr, strand, site1["pos"], site2["pos"])
                pair_cDNA_sum = site1["cDNA_sum"] + site2["cDNA_sum"]
                site_pairs.append((pair_cDNA_sum, pair_type, site1, site2))
                site_pairs_stats[pair_type]+=1

            site_pairs.sort(reverse=True)

            for (cDNA_pair, pair_type, site1, site2) in site_pairs:

                if pair_type in considered_types:
                    continue

                # todo: check if this makes sense?
                major = max(site1["cDNA_sum"], site2["cDNA_sum"])
                minor = min(site1["cDNA_sum"], site2["cDNA_sum"])
                if minor < major*0.1:
                    continue

                considered_types.add(pair_type)

                #site1 = L[0][1]
                #site2 = L[1][1]
                #pair_type = apa.polya.annotate_pair(comps.species, chr, strand, site1["pos"], site2["pos"])

                # how many tandem sites are we loosing?
                #if len(site_pairs)>1:
                #    if "tandem" in site_pairs_stats and pair_type!="tandem":
                #        print site_pairs_stats, pair_type # all sites types vs. major pair type


                if (site1["pos"]<site2["pos"] and strand=="+") or (site1["pos"]>site2["pos"] and strand=="-"):
                    up_site = site1
                    down_site = site2
                else:
                    up_site = site2
                    down_site = site1

                up_control = []
                up_test = []
                for (rshort, _) in replicates:
                    if rshort.startswith("c"):
                        up_control.append(up_site[rshort])
                    else:
                        up_test.append(up_site[rshort])

                down_control = []
                down_test = []
                for (rshort, _) in replicates:
                    if rshort.startswith("c"):
                        down_control.append(down_site[rshort])
                    else:
                        down_test.append(down_site[rshort])

                up_seq = pybio.genomes.seq(comps.species, chr, strand, site1["pos"], start=-60, stop=100)
                down_seq = pybio.genomes.seq(comps.species, chr, strand, site2["pos"], start=-60, stop=100)

                _, up_vector = pybio.sequence.search(up_seq, ["TGT", "GTG"])
                up_vector = pybio.sequence.filter(up_vector, hw=25, hwt=17)
                _, down_vector = pybio.sequence.search(down_seq, ["TGT", "GTG"])
                down_vector = pybio.sequence.filter(down_vector, hw=25, hwt=17)

                row = [chr, strand, gene_locus, gene_id, gene["gene_name"], gene["gene_biotype"], len(sites)]
                row.append(up_site["pos"])
                row.append(sum(up_test+up_control))
                row.append(sum(up_vector))
                row.append(down_site["pos"])
                row.append(sum(down_test+down_control))
                row.append(sum(down_vector))

                row.append(";".join(str(x) for x in up_control))
                row.append(sum(up_control))
                row.append(";".join(str(x) for x in down_control))
                row.append(sum(down_control))
                row.append(";".join(str(x) for x in up_test))
                row.append(sum(up_test))
                row.append(";".join(str(x) for x in down_test))
                row.append(sum(down_test))

                try:
                    pc = float(sum(up_control))/sum(up_control + down_control) - float(sum(up_test))/sum(up_test + down_test)
                except:
                    pc = 0
                row.append("%.5f" % pc)

                f = fisher.pvalue(sum(up_control), sum(up_test), sum(down_control), sum(down_test))
                pvalue = f.two_tail
                row.append("%.5f" % pvalue)
                pair_type = apa.polya.annotate_pair(comps.species, chr, strand, up_site["pos"], down_site["pos"])
                row.append(pair_type)
                results.append(row)

    results = sorted(results, key=lambda x: abs(float(x[-3])), reverse=True)
    results = sorted(results, key=lambda x: float(x[-2]))
    results = sorted(results, key=lambda x: x[-1], reverse=True)
    for row in results:
        f_pairs.write("\t".join([str(x) for x in row]) + "\n")
    f_pairs.close()


def distance_hist(comps_id):
    pairs_filename = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
    hist = []
    f = open(pairs_filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        proximal = int(data["siteup_pos"])
        distal = int(data["sitedown_pos"])
        d = abs(proximal-distal)
        if d<10000:
            hist.append(d)
        r = f.readline()
    f.close()

    print hist

    import numpy as np
    import pylab as P
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tkr

    def func(x, pos):  # formatter function takes tick label and tick position
       s = '{:0,d}'.format(int(x))
       return s

    n, bins, patches = P.hist(hist, 20, histtype='stepfilled')
    P.setp(patches, 'facecolor', 'gray', 'alpha', 0.75)
    fig_fname = os.path.join(apa.path.comps_folder, comps_id, "%s.proximal_distal.png" % comps_id)
    P.title(comps_id)
    P.ylim(bottom=0)
    P.xlabel("proximal to distal [nt]")
    P.ylabel("number of genes")

    a_format = tkr.FuncFormatter(func)
    axes = plt.gca()
    axes.xaxis.set_major_formatter(a_format)
    axes.xaxis.set_major_formatter(a_format)
    axes.spines['bottom'].set_alpha(0.5)
    axes.spines['top'].set_alpha(0.5)
    axes.spines['right'].set_alpha(0.5)
    axes.spines['left'].set_alpha(0.5)
    P.savefig(fig_fname)


def utr_boxplot(study_id, comps_id):
    comps = read_comps(study_id, comps_id) # study data

    replicates = []
    expression = {} # keys = c1, c2, c3, t1, t2, t3...items = bedgraph files

    for (comp_id, experiments, comp_name) in comps.test+comps.control:
        key = "%s:%s" % (comp_id, comp_name)
        print "reading: ", key
        expression[comp_id] = pybio.data.Bedgraph()
        replicates.append((comp_id, key)) # (short_name, long_name)
        for id in experiments:
            lib_id = id[:id.rfind("_")]
            exp_id = int(id.split("_")[-1][1:])
            e_filename = apa.path.e_filename(lib_id, exp_id)
            expression[comp_id].load(e_filename)
        print

    replicates.sort()

    # read study gene sites
    genes = {}
    expression_genes_pair_tandem_filename = apa.path.comps_expression_filename(study_id, comps_id, filetype="genes.pair")
    f = open(expression_genes_pair_tandem_filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        genes[data["gene_id"]] = data
        r = f.readline()
    f.close()

    sitec_box1 = [] # control box 1
    sitec_box2 = [] # control box 2
    sitet_box1 = [] # test box 1
    sitet_box2 = [] # test box 2
    for gid, gene_data in genes.items():
        gene = apa.polya.get_gene(comps.species, gid)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)
        row2 = []
        cDNA_total = 0
        if gene_data.get("siteup_pos", None)==None:
            # no data for this gene
            continue
        site_positions = [int(gene_data["siteup_pos"]), int(gene_data["sitedown_pos"])]

        if gene_data["pair_type"]!="tandem":
            continue

        site1c = 0
        site2c = 0
        site1t = 0
        site2t = 0
        for (comp_id, _, _) in comps.control:
            bg = expression[comp_id]
            site1c += bg.get_value(chr, strand, site_positions[0], db="raw")
            site2c += bg.get_value(chr, strand, site_positions[1], db="raw")
        for (comp_id, _, _) in comps.test:
            bg = expression[comp_id]
            site1t += bg.get_value(chr, strand, site_positions[0], db="raw")
            site2t += bg.get_value(chr, strand, site_positions[1], db="raw")

        # normalize by number of control and test
        #site1c /= len(comps.control)
        #site2c /= len(comps.control)
        #site1t /= len(comps.test)
        #site2t /= len(comps.test)

        # sum
        sitec_all = site1c + site2c
        sitet_all = site1t + site2t

        # box
        if sitec_all>50 and sitet_all>50:
            sitec_box1.append(site1c/float(sitec_all))
            sitec_box2.append(site2c/float(sitec_all))
            sitet_box1.append(site1t/float(sitet_all))
            sitet_box2.append(site2t/float(sitet_all))

    matplotlib.rcParams['axes.labelsize'] = 12
    matplotlib.rcParams['axes.titlesize'] = 12
    matplotlib.rcParams['xtick.labelsize'] = 12
    matplotlib.rcParams['ytick.labelsize'] = 12
    matplotlib.rcParams['legend.fontsize'] = 10
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.3
    matplotlib.rcParams['legend.frameon'] = 'False'

    x = numpy.array([1000, 2000])
    datac = [sitec_box1, sitec_box2]
    datat = [sitet_box1, sitet_box2]
    meansc = [numpy.mean(el) for el in datac]
    meanst = [numpy.mean(el) for el in datat]
    xlabels = ["site.up", "site.down"]

    plt.figure(figsize=(10, 5))
    a = plt.axes([0.1, 0.1, 0.85, 0.8])
    a.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    bpc = plt.boxplot(datac, 1, '', positions=x-100, widths=100)
    plt.setp(bpc['boxes'], color='black')
    plt.setp(bpc['whiskers'], color='black')
    plt.setp(bpc['fliers'], color='red', marker='+')
    for line in bpc['boxes']:
        line.set_linewidth(0.6)
    for line in bpc['whiskers']:
        line.set_linewidth(0.5)

    bpt = plt.boxplot(datat, 1, '', positions=x+100, widths=100)
    plt.setp(bpt['boxes'], color='blue')
    plt.setp(bpt['whiskers'], color='blue')
    plt.setp(bpt['fliers'], color='red', marker='+')
    for line in bpt['boxes']:
        line.set_linewidth(0.6)
    for line in bpt['whiskers']:
        line.set_linewidth(0.5)
    plt.xlim(500, 2500)
    plt.xticks([1000, 2000], xlabels)
    plt.savefig(os.path.join(apa.path.study_folder, study_id, comps_id, "%s.%s" % (comps_id, "tandem.utr.svg")))

def make_fasta(comps_id):
    comps = apa.comps.read_comps(comps_id)

    stats_reg = {}
    for reg in ["siteup.c", "siteup.e", "siteup.r", "sitedown.c", "sitedown.e", "sitedown.r"]:
        stats_reg[reg] = 0

    # fasta
    fasta_folder = os.path.join(apa.path.comps_folder, comps_id, "fasta")
    if os.path.exists(fasta_folder):
        shutil.rmtree(fasta_folder)
    os.makedirs(fasta_folder)

    par = {0:0.1, 1:0.05, 2:2, 3:"tandem"}
    genes = {}
    input_file = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
    f = open(input_file, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        genes[data["gene_id"]] = data
        r = f.readline()
    f.close()

    for gene_id, data in genes.items():
        siteup_pos = int(data["siteup_pos"])
        sitedown_pos = int(data["sitedown_pos"])
        pc = float(data["pc"])
        fisher = float(data["fisher"])
        pair_type = data["pair_type"]
        if par[3]!=pair_type:
            continue
        site_distance = abs(siteup_pos-sitedown_pos)+1
        if par[2]==0 and site_distance>200:
            continue
        if par[2]==1 and (site_distance<200 or site_distance>450):
            continue
        if par[2]==2 and site_distance<450:
            continue
        gene = apa.polya.get_gene(comps.species, gene_id)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]

        reg = None
        if pc>=0 and abs(pc)>par[0] and fisher<par[1]:
            reg = "r"
        if pc<0 and abs(pc)>par[0] and fisher<par[1]:
            reg = "e"
        if fisher>0.5:
            reg = "c"
        if reg==None:
            continue

        siteup_reg = reg
        sitedown_reg = {"e":"r", "r":"e", "c":"c"}[reg]

        for (site_reg, site_type, site_pos) in [(siteup_reg, "siteup", siteup_pos), (sitedown_reg, "sitedown", sitedown_pos)]:
            reg_key = "%s.%s" % (site_type, site_reg)
            stats_reg[reg_key] = stats_reg[reg_key]+1
            seq = pybio.genomes.seq(comps.species, chr, strand, site_pos, start=-50, stop=50)
            fasta_filename = os.path.join(apa.path.comps_folder, comps_id, "fasta", "%s.fasta" % (reg_key))
            if not os.path.exists(fasta_filename):
                f = open(fasta_filename, "wt")
            else:
                f = open(fasta_filename, "at")
            f.write(">%s.%s\n%s\n" % (reg_key, stats_reg[reg_key], seq))
            f.close()
