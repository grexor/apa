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
import copy

minor_major_thr = 0.05 # for a site to be considered, it's expression needs to be >5% of the max site inside the gene
clip_interval = (-100, 100) # interval around sites for binding data

class Comps:
    def __init__(self, comps_id=None):
        self.comps_id = comps_id
        self.test = []
        self.control = []
        self.test_name = ""
        self.control_name = ""
        self.species = ""
        self.CLIP = []
        self.cDNA_thr = 5 # at least cDNA
        self.presence_thr = 2.0 # for at least half of experiments
        self.pc_thr = 0.1
        self.fisher_thr = 0.1
        self.control_thr = 0.025
        self.pair_dist = 450
        self.exp_data = {}
        self.polya_db = ""
        self.poly_type = ["strong", "weak"] # strong, weak, less, noclass
        self.site_selection = "CLIP" # CLIP / APA
        self.deepbind = None
        self.rnamaps = []
        self.use_FDR = False
        self.ignore_genes = []
        self.exclusive_genes = []
        self.db_type="cs" # cs = cleavage site, pas = polyadenylation signal, cs = default

    def __str__(self):
        print "comps_id = %s" % (self.comps_id)
        for test_id, test_data in self.test.items():
            print "%s %s" % (test_id, test_data)
        for control_id, control_data in self.control.items():
            print "%s %s" % (control_id, control_data)
        print "CLIP = %s" % self.CLIP
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
        if r[0].startswith("use_FDR:"):
            if r[0].split("use_FDR:")[1] in ["True", "true", "yes", "1"]:
                comps.use_FDR = True
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
            fname = r[0].split("iCLIP:")[1]
            comps.CLIP.append(fname)
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
        if r[0].startswith("site_selection"):
            comps.site_selection = r[0].split("site_selection:")[1]
            r = f.readline()
            continue
        if r[0].startswith("polya_db:"):
            comps.polya_db = r[0].split("polya_db:")[1]
            r = f.readline()
            continue
        if r[0].startswith("poly_type:"):
            comps.poly_type = eval(r[0].split("poly_type:")[1])
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
        if r[0].startswith("ignore_genes:"):
            comps.ignore_genes = r[0].split("ignore_genes:")[1].split(",")
            r = f.readline()
            continue
        if r[0].startswith("exclusive_genes:"):
            comps.exclusive_genes = eval(r[0].split("exclusive_genes:")[1])
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
    for (_, exp, _) in comps.test+comps.control:
        if type(exp)==list:
            for exp_id in exp:
                all_exp.update(exp)
        else:
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
    if len(comps.CLIP)>0:
        f.write("CLIP:%s" % str(comps.CLIP))
    f.close()
    return True

def process_comps(comps_id, map_id=1):

    if not os.path.exists(os.path.join(apa.path.comps_folder, comps_id, "%s.config" % comps_id)):
        print "%s.config missing, quitting" % (comps_id)
        sys.exit(1)

    # clean
    assert(len(apa.path.comps_folder)>0)
    files = glob.glob(os.path.join(apa.path.comps_folder, comps_id, "*"))
    for f in files:
        if not f.endswith(".config") and not f.startswith("docs"):
            print "removing: %s" % f
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)

    comps = read_comps(comps_id) # study data
    pybio.genomes.load(comps.species)

    # load CLIP data if available
    clip = {}
    for clip_name in comps.CLIP:
        clip[clip_name] = pybio.data.Bedgraph(os.path.join(apa.path.iCLIP_folder, clip_name))

    # if there is a polya-db specified in the comparison, load the positions into the filter
    # (strong, weak, less)
    poly_filter = {}
    polydb = pybio.data.Bedgraph()
    for poly_type in comps.poly_type:
        polydb.load(apa.path.polyadb_filename(comps.polya_db, poly_type=poly_type, filetype="bed"))

    replicates = []
    expression = {} # keys = c1, c2, c3, t1, t2, t3...items = bedgraph files

    for (comp_id, experiments, comp_name) in comps.control+comps.test:
        key = "%s:%s" % (comp_id, comp_name)
        expression[comp_id] = pybio.data.Bedgraph()
        replicates.append((comp_id, key)) # (short_name, long_name)
        for id in experiments:
            print "%s: +%s" % (comp_id, id)
            lib_id = id[:id.rfind("_")]
            exp_id = int(id.split("_")[-1][1:])
            e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=comps.polya_db)
            expression[comp_id].load(e_filename)
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
                    e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=comps.polya_db)
                elif comps.db_type=="pas":
                    e_filename = apa.path.e_filename(lib_id, exp_id, filetype="pas", map_id=map_id, poly_id=comps.polya_db)
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
                    # filter poly-A positions by type (strong, weak, etc)
                    if polydb.get_value(chr, strand, pos)!=0:
                        valid_positions.add(pos)
                positions[chr][strand] = positions[chr][strand].union(valid_positions)

    # organize polya sites / pas signals inside genes
    gsites = {}

    num_sites = 0
    num_sites_genes = 0
    num_genes = set()
    num_sites_genes_expressed = 0
    num_sites_genes_expressed_final = 0

    for chr, strand_data in positions.items():
        for strand, pos_set in strand_data.items():
            pos_set = list(pos_set)
            for pos in pos_set:
                num_sites += 1
                gene_up, gene_id, gene_down, gene_interval = apa.polya.annotate_position(comps.species, chr, strand, pos)
                if gene_id==None: # only consider polya sites inside genes
                    continue
                num_sites_genes += 1
                num_genes.add(gene_id)
                sites = gsites.get(gene_id, {})
                expression_vector = []
                site_data = {"pos":pos, "gene_interval":list(gene_interval)} # store position and gene_interval (start, stop, exon/intron), clip binding

                # get clip data
                for clip_name in comps.CLIP:
                    site_data[clip_name] = clip[clip_name].get_region("chr"+chr, strand, pos, start=clip_interval[0], stop=clip_interval[1])

                control_sum = 0
                test_sum = 0
                cDNA_sum = 0
                for (rshort, _) in replicates:
                    bg = expression[rshort]
                    cDNA = bg.get_value(chr, strand, pos)
                    cDNA_sum += cDNA
                    expression_vector.append(cDNA)
                    if rshort.startswith("c"):
                        control_sum += cDNA
                    if rshort.startswith("t"):
                        test_sum += cDNA
                    site_data[rshort] = cDNA

                # OLD: filter lowly expressed sites
                #expression_vector = [1 if cDNA>=comps.cDNA_thr else 0 for cDNA in expression_vector]
                #if sum(expression_vector) < len(expression_vector)/comps.presence_thr:
                #    continue

                # filter lowly expressed sites
                if (control_sum<10) and (test_sum<10):
                    continue

                site_data["cDNA_sum"] = int(cDNA_sum)
                sites[pos] = site_data
                gsites[gene_id] = sites
                num_sites_genes_expressed += 1

    # filter out sites that have expression < minor_major_thr of maximally expressed site
    for gene_id, sites in gsites.items():
        max_exp = 0
        for pos, site_data in sites.items():
            max_exp = max(max_exp, site_data["cDNA_sum"])
        for pos, site_data in list(sites.items()): # python3 safe
            if site_data["cDNA_sum"] < (minor_major_thr * max_exp):
                del sites[pos]

    num_sites_per_gene = {}
    for gene_id, sites in gsites.items():
        num_sites_per_gene[len(sites.keys())] = num_sites_per_gene.get(len(sites.keys()), 0) + 1
        num_sites_genes_expressed_final += len(sites.keys())

    f_stats = open(os.path.join(apa.path.comps_folder, comps_id, "%s_annotation_stats.txt" % comps_id), "wt")
    f_stats.write("-----------site stats---------------\n")
    f_stats.write("%s polyA sites\n" % num_sites)
    f_stats.write("    %s polyA sites assigned to %s genes\n" % (num_sites_genes, len(num_genes)))
    f_stats.write("        %s polyA sites expressed enough across experiments\n" % num_sites_genes_expressed)
    f_stats.write("            %s final polyA sites kept that had >5%% expression on gene level\n" % num_sites_genes_expressed_final)
    f_stats.write("-----------site stats (end)---------\n")
    f_stats.write("\n")
    f_stats.write("-----------site per gene -----------\n")
    sum_temp_sites = 0
    sum_temp_genes = 0
    for s_gene, v in num_sites_per_gene.items():
        f_stats.write("%s genes with %s site(s)\n" % (v, s_gene))
        sum_temp_sites += s_gene*v
        sum_temp_genes += v
    f_stats.write("total of %s polyA sites assigned to %s genes\n" % (sum_temp_sites, sum_temp_genes))
    f_stats.write("-----------site per gene (end)------\n")
    f_stats.close()

    # draw num sites per gene
    # =======================
    import matplotlib
    matplotlib.use("Agg", warn=False)
    import matplotlib.pyplot as plt
    import math
    import gzip
    from matplotlib import cm as CM
    import matplotlib.patches as mpatches
    import matplotlib.ticker as mticker
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.colors as mcolors
    matplotlib.rcParams['axes.labelsize'] = 14
    matplotlib.rcParams['axes.titlesize'] = 14
    matplotlib.rcParams['xtick.labelsize'] = 14
    matplotlib.rcParams['ytick.labelsize'] = 14
    matplotlib.rcParams['legend.fontsize'] = 14
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.4
    matplotlib.rcParams['legend.frameon'] = 'False'
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., height+100, format(int(height), ','), ha='center', va='bottom', fontdict={"size":14})
    fig, ax = plt.subplots(1, 1, figsize=(12, 3))
    total_sites, total_genes = 0, 0
    x = [1,2,3,4,5,6,7,8,9,10]
    y = {}
    for e in x:
        y[e] = 0
    for x1 in sorted(num_sites_per_gene.keys()):
        if x1<10:
            y[x1] += num_sites_per_gene[x1]
        else:
            y[10] += num_sites_per_gene[x1]
        total_sites += num_sites_per_gene[x1] * x1
        total_genes += num_sites_per_gene[x1]
    y = [y[e] for e in x]
    bar1 = ax.bar(x, y, 0.5, label='number of polyA sites', color='lightgray')
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    plt.xlim(1-0.25, 10+0.5)
    plt.xlabel('number of polyA sites'); plt.ylabel('number of genes')
    plt.xticks([e+0.25 for e in x], [1,2,3,4,5,6,7,8,9,">=10"])
    autolabel(bar1)
    plt.title("%s polyA sites annotated to %s genes" % (format(int(total_sites), ','), format(int(total_genes), ',')))
    plt.tight_layout()
    plt.savefig(os.path.join(apa.path.comps_folder, comps_id, "%s_sites_per_gene.pdf" % comps_id))
    plt.savefig(os.path.join(apa.path.comps_folder, comps_id, "%s_sites_per_gene.png" % comps_id), transparent=True)
    plt.close()
    # =======================

    # expression gene level
    fname = apa.path.comps_expression_filename(comps_id)
    f_genes = open(fname, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype"]
    header += ["sites_position_cDNA", "sites_num", "cDNA_sum"]
    # would be nice to have entire gene clip binding, but takes some time
    #for clip_name in comps.CLIP:
    #    header.append("clip:%s" % clip_name)
    for (rshort, rlong) in replicates:
        header.append(rlong)
    f_genes.write("\t".join(header)+"\n")
    for gene_id, sites in gsites.items():
        gene = apa.polya.get_gene(comps.species, gene_id)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_len = gene_stop-gene_start+1
        assert(gene_len>0)
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
        # would be nice to have entire gene clip binding, but takes some time
        #for clip_name in comps.CLIP:
        #    row.append(clip[clip_name].get_region("chr"+chr, strand, gene_start, stop=gene_len))
        row += row_comp
        f_genes.write("\t".join([str(x) for x in row]) + "\n")
    f_genes.close()

    # expression site level
    gene_sites = {}
    fname = apa.path.comps_expression_filename(comps_id, filetype="sites")
    f_sites = open(fname, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "site_pos", "gene_interval"]
    for clip_name in comps.CLIP:
        header.append("clip:%s" % clip_name)
    header.append("cDNA_sum")
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
        for site_pos, site_data in sites.items():
            row_1 = [chr, strand, gene_locus, gene_id, gene["gene_name"], gene["gene_biotype"], site_pos, site_data["gene_interval"]]
            for clip_name in comps.CLIP:
                row_1.append(site_data[clip_name])
            row_2 = []
            cDNA_sum = 0
            for (rshort, rlong) in replicates:
                val = site_data.get(rshort, 0)
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
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "num_sites", "proximal_pos", "proximal_exp", "proximal_UG", "distal_pos", "distal_exp", "distal_UG", "s1", "s2"]
    header.append("proximal_control")
    header.append("proximal_control_sum")
    header.append("distal_control")
    header.append("distal_control_sum")
    header.append("proximal.test")
    header.append("proximal_test_sum")
    header.append("distal.test")
    header.append("distal_test_sum")

    header += ["pc", "fisher", "pair_type"]

    bg_selected_sites_control = pybio.data.Bedgraph()
    bg_selected_sites_test = pybio.data.Bedgraph()

    num_genes = {}
    f_pairs.write("\t".join(header)+"\n")
    results = []
    for gene_id, sites in gsites.items():
        gene = apa.polya.get_gene(comps.species, gene_id)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)

        already_considered_types = set()

        if len(sites)<2:
            continue

        if len(comps.CLIP)>0 and comps.site_selection=="CLIP":
            L = [(sites[pos][comps.CLIP[0]], sites[pos]["cDNA_sum"], sites[pos]) for pos in sites.keys()] # we take the first CLIP file and use it to determine regulated sites
        else:
            L = [(0, sites[pos]["cDNA_sum"], sites[pos]) for pos in sites.keys()]

        # determine major / minor sites
        S_clip = copy.deepcopy(L)
        S_clip.sort(key=lambda x: x[0], reverse=True) # sort by clip_binding
        S_exp = copy.deepcopy(L)
        S_exp.sort(key=lambda x: x[1], reverse=True) # sort by expression

        if S_clip[0][0]>0:
            major = S_clip[0][-1] # take position of most bound site
            minor = S_exp[0][-1] # take position of most expressed site
            if major["pos"]==minor["pos"]:
                minor = S_exp[1][-1]
        else:
            major = S_exp[0][-1] # most expressed
            minor = S_exp[1][-1] # second most expressed

        pair_type = apa.polya.annotate_pair(comps.species, chr, strand, major["pos"], minor["pos"]) # what kind of pair is it?
        num_genes[pair_type] = num_genes.get(pair_type, 0) + 1
        minor_sites = []
        if pair_type=="tandem":
            for site_pos, site_data in sites.items():
                if site_data["gene_interval"]==major["gene_interval"] and major["pos"]!=site_data["pos"]:
                    minor_sites.append(site_data)
            assert(len(minor_sites)>0) # tandem? at least 1 other site in same exon
        else:
            for site_pos, site_data in sites.items():
                if major["pos"]!=site_data["pos"]:
                    minor_sites.append(site_data)
            assert(len(minor_sites)>0) # skipped or composite? at least 1 other site in gene

        # sum up counts at minor sites
        # first clear the minor site values
        minor["cDNA_sum"] = 0
        for (rshort, rlong) in replicates:
            minor[rshort] = 0
        # then add the values from minor_sites list
        for (rshort, rlong) in replicates:
            for site_data in minor_sites:
                minor[rshort] += site_data[rshort]
        # finally sum up cDNA
        cDNA_sum = 0
        for (rshort, rlong) in replicates:
            cDNA_sum += minor[rshort]
        minor["cDNA_sum"] = cDNA_sum

        if (major["pos"]<minor["pos"] and strand=="+") or (major["pos"]>minor["pos"] and strand=="-"):
            proximal_site = major
            distal_site = minor
        else:
            proximal_site = minor
            distal_site = major

        s1, s2 = get_s1_s2(gene_id, chr, strand, comps.species, proximal_site["pos"], distal_site["pos"], pair_type)

        proximal_control = []
        proximal_test = []
        for (rshort, _) in replicates:
            if rshort.startswith("c"):
                proximal_control.append(proximal_site[rshort])
            else:
                proximal_test.append(proximal_site[rshort])

        distal_control = []
        distal_test = []
        for (rshort, _) in replicates:
            if rshort.startswith("c"):
                distal_control.append(distal_site[rshort])
            else:
                distal_test.append(distal_site[rshort])

        # update bedGraph for selected sites
        bg_selected_sites_control.set_value("chr"+chr, strand, proximal_site["pos"], sum(proximal_control))
        bg_selected_sites_test.set_value("chr"+chr, strand, proximal_site["pos"], sum(proximal_test))
        bg_selected_sites_control.set_value("chr"+chr, strand, distal_site["pos"], sum(distal_control))
        bg_selected_sites_test.set_value("chr"+chr, strand, distal_site["pos"], sum(distal_test))

        proximal_seq = pybio.genomes.seq(comps.species, chr, strand, proximal_site["pos"], start=-60, stop=100)
        distal_seq = pybio.genomes.seq(comps.species, chr, strand, distal_site["pos"], start=-60, stop=100)
        _, proximal_vector = pybio.sequence.search(proximal_seq, ["TGT", "GTG"])
        proximal_vector = pybio.sequence.filter(proximal_vector, hw=25, hwt=17)
        _, distal_vector = pybio.sequence.search(distal_seq, ["TGT", "GTG"])
        distal_vector = pybio.sequence.filter(distal_vector, hw=25, hwt=17)

        row = [chr, strand, gene_locus, gene_id, gene["gene_name"], gene["gene_biotype"], len(sites)]
        row.append(proximal_site["pos"])
        row.append(sum(proximal_test+proximal_control))
        row.append(sum(proximal_vector))
        row.append(distal_site["pos"])
        row.append(sum(distal_test+distal_control))
        row.append(sum(distal_vector))

        row.append(s1)
        row.append(s2)

        row.append(";".join(str(x) for x in proximal_control))
        row.append(sum(proximal_control))
        row.append(";".join(str(x) for x in distal_control))
        row.append(sum(distal_control))
        row.append(";".join(str(x) for x in proximal_test))
        row.append(sum(proximal_test))
        row.append(";".join(str(x) for x in distal_test))
        row.append(sum(distal_test))

        try:
            pc = float(sum(proximal_control))/sum(proximal_control + distal_control) - float(sum(proximal_test))/sum(proximal_test + distal_test)
        except:
            pc = 0
        row.append("%.5f" % pc)

        f = fisher.pvalue(sum(proximal_control), sum(proximal_test), sum(distal_control), sum(distal_test))
        pvalue = f.two_tail
        row.append("%.5f" % pvalue)
        pair_type = apa.polya.annotate_pair(comps.species, chr, strand, proximal_site["pos"], distal_site["pos"])
        row.append(pair_type)
        results.append(row)

    results = sorted(results, key=lambda x: abs(float(x[-3])), reverse=True)
    results = sorted(results, key=lambda x: float(x[-2]))
    results = sorted(results, key=lambda x: x[-1], reverse=True)
    for row in results:
        f_pairs.write("\t".join([str(x) for x in row]) + "\n")
    f_pairs.close()

    # write expression_proximal file
    f = open(os.path.join(apa.path.comps_folder, comps_id, "%s.expression_proximal.tab" % comps_id), "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "num_sites", "proximal_pos", "proximal_exp"]
    i = 0
    for (rshort, _) in replicates:
        if rshort.startswith("c"):
            i += 1
            header.append("c%s" % i)
    i = 0
    for (rshort, _) in replicates:
        if rshort.startswith("t"):
            i += 1
            header.append("t%s" % i)
    header.append("csum")
    header.append("tsum")
    header.append("pair_type")
    f.write("\t".join(header)+"\n")
    for row in results:
        r = row[:9]
        r = r + row[15].split(";")
        r = r + row[19].split(";")
        r.append(row[16])
        r.append(row[20])
        r.append(row[-1])
        f.write("\t".join([str(x) for x in r]) + "\n")
    f.close()

    # write expression_distal file
    f = open(os.path.join(apa.path.comps_folder, comps_id, "%s.expression_distal.tab" % comps_id), "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "num_sites", "distal_pos", "distal_exp"]
    i = 0
    for (rshort, _) in replicates:
        if rshort.startswith("c"):
            i += 1
            header.append("c%s" % i)
    i = 0
    for (rshort, _) in replicates:
        if rshort.startswith("t"):
            i += 1
            header.append("t%s" % i)
    header.append("csum")
    header.append("tsum")
    header.append("pair_type")
    f.write("\t".join(header)+"\n")
    for row in results:
        r = row[:7]
        r.append(row[10])
        r.append(row[11])
        r = r + row[17].split(";")
        r = r + row[21].split(";")
        r.append(row[18])
        r.append(row[22])
        r.append(row[-1])
        f.write("\t".join([str(x) for x in r]) + "\n")
    f.close()

    if comps.use_FDR:
        pybio.utils.FDR_tab(pairs_filename, "fisher")

    print "----------stats pairs------------"
    print num_genes.items()
    print "----------stats pairs (end)------"

    # save selected sites bedGraphs
    bg_selected_sites_control_fname = os.path.join(beds_folder, "%s_control_selected.bed" % comps_id)
    bg_selected_sites_test_fname = os.path.join(beds_folder, "%s_test_selected.bed" % comps_id)

    bg_selected_sites_control.save(bg_selected_sites_control_fname, track_id="%s_control_selected" % comps_id)
    bg_selected_sites_test.save(bg_selected_sites_test_fname, track_id="%s_test_selected" % comps_id)

def get_s1_s2(gene_id, chr, strand, genome, proximal_pos, distal_pos, pair_type):
    pybio.genomes.load(genome)
    gene = apa.polya.get_gene(genome, gene_id)
    _, _, _, gene_interval = apa.polya.annotate_position(genome, chr, strand, proximal_pos)
    proximal_site = {"pos":proximal_pos, "gene_interval":list(gene_interval)}
    _, _, _, gene_interval = apa.polya.annotate_position(genome, chr, strand, distal_pos)
    distal_site = {"pos":distal_pos, "gene_interval":list(gene_interval)}

    s1 = None
    s2 = None
    # https://docs.google.com/drawings/d/1_m4iZ1c9YwKI-NOWMSCSdg2IzEGedj-NaMTIlHqirc0/edit

    if pair_type=="tandem":
        interval = proximal_site["gene_interval"]
        intervals = gene["gene_intervals"]
        interval_index = intervals.index(interval)
        interval_upstream_index = (interval_index - 1) if strand=="+" else (interval_index + 1)
        if interval_upstream_index==-1 or interval_upstream_index>=len(intervals):
            interval_upstream = None
        else:
            interval_upstream = intervals[interval_upstream_index]
            # some genes are split, e.g. ENSG00000163684 (because of shorter gene priority when encountering overlapping genes)
            # this causes the gene not to have a structure of oioioi; a situation like oioo can happen
            # do not consider those cases
            if interval_upstream[-1]!="i":
                interval_upstream = None
        if interval_upstream!=None:
            if strand=="+":
                s1 = interval_upstream[0]
                s2 = interval_upstream[1]
            else:
                s1 = interval_upstream[1]
                s2 = interval_upstream[0]

    if pair_type=="composite":
        interval = proximal_site["gene_interval"]
        assert(interval[-1]=="i")
        if strand=="+":
            s1 = interval[0]
            s2 = interval[1]
        else:
            s1 = interval[1]
            s2 = interval[0]

    if pair_type=="skipped":
        if strand=="+":
            s1 = proximal_site["gene_interval"][0]
            s2 = distal_site["gene_interval"][0]
        else:
            s1 = proximal_site["gene_interval"][1]
            s2 = distal_site["gene_interval"][1]

    return s1, s2

def distance_hist(comps_id):
    pairs_filename = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
    hist = []
    f = open(pairs_filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        proximal = int(data["proximal_pos"])
        distal = int(data["distal_pos"])
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
            e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=comps.polya_db)
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
        if gene_data.get("proximal_pos", None)==None: # no data for this gene
            continue
        site_positions = [int(gene_data["proximal_pos"]), int(gene_data["distal_pos"])]

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
    for reg in ["proximal.c", "proximal.e", "proximal.r", "distal.c", "distal.e", "distal.r"]:
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
        proximal_pos = int(data["proximal_pos"])
        distal_pos = int(data["distal_pos"])
        pc = float(data["pc"])
        fisher = float(data["fisher"])
        pair_type = data["pair_type"]
        if par[3]!=pair_type:
            continue
        site_distance = abs(proximal_pos-distal_pos)+1
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

        proximal_reg = reg
        distal_reg = {"e":"r", "r":"e", "c":"c"}[reg]

        for (site_reg, site_type, site_pos) in [(proximal_reg, "proximal", proximal_pos), (distal_reg, "distal", distal_pos)]:
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
