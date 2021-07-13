import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import math
from matplotlib import cm as CM
import numpy
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import pybio
import apa
import os
import shutil
import json
import glob
import sys
import operator
import itertools
from collections import Counter
import copy

minor_major_thr = 0.05 # for a site to be considered, it's expression needs to be >5% of the max site inside the gene

class Comps:

    def __init__(self, comps_id=None):
        self.comps_id = comps_id
        self.test = []
        self.control = []
        self.test_name = ""
        self.control_name = ""
        self.species = ""
        self.CLIP = []
        self.name = ""
        self.notes = ""
        self.access = ""
        self.authors = ""
        self.owner = []
        self.cDNA_thr = 5 # at least cDNA
        self.presence_thr = 2.0 # for at least half of experiments
        self.significance_thr = 0.05
        self.pair_dist = 0
        self.exp_data = {}
        self.polya_db = ""
        self.poly_type = ["strong", "weak"] # strong, weak, less, noclass
        self.site_selection = "DEXSEQ" # possibility to do site selection some other way in the future
        self.choose_function = "sum" # how are sites evaluated from iCLIP data, by "sum" in the clip_interval or by "closest" in the "clip_interval"
        self.deepbind = None
        self.cluster_image_w = 1000
        self.cluster_image_h = 1000
        self.clip_interval = (-100, 100) # interval around sites for binding data
        self.rnamaps = []
        self.ignore_genes = []
        self.exclusive_genes = []
        self.genome = ""
        self.method = ""
        self.db_type="cs" # cs = cleavage site, pas = polyadenylation signal, cs = default
        self.analysis_type = "apa"
        self.status = "complete" # processing, complete, default is complete since we dont have this parameters for all analysis
        if comps_id!=None:
            self.read_comps(comps_id)

    def __str__(self):
        print("comps_id = %s" % (self.comps_id))
        for test_id, test_data in self.test.items():
            print("%s %s" % (test_id, test_data))
        for control_id, control_data in self.control.items():
            print("%s %s" % (control_id, control_data))
        print("CLIP = %s" % self.CLIP)
        return ""

    def read_comps(self, comps_id):
        config_file = apa.path.comps_config_filename(comps_id)
        if not os.path.exists(config_file):
            return
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
            if r[0].startswith("genome"):
                self.genome = r[0].split("genome:")[1]
                r = f.readline()
                continue
            if r[0].startswith("analysis_type"):
                self.analysis_type = r[0].split("analysis_type:")[1]
                r = f.readline()
                continue
            if r[0].startswith("method"):
                self.method = r[0].split("method:")[1]
                r = f.readline()
                continue
            if r[0].startswith("status"):
                self.status = r[0].split("status:")[1]
                r = f.readline()
                continue
            if r[0].startswith("cluster_image_w"):
                self.cluster_image_w = int(r[0].split("cluster_image_w:")[1])
                r = f.readline()
                continue
            if r[0].startswith("cluster_image_h"):
                self.cluster_image_h = int(r[0].split("cluster_image_h:")[1])
                r = f.readline()
                continue
            if r[0].startswith("access:"):
                self.access = str(r[0].split("access:")[1]).split(",")
                r = f.readline()
                continue
            if r[0].startswith("owner:"):
                self.owner = str(r[0].split("owner:")[1]).split(",")
                r = f.readline()
                continue
            if r[0].startswith("name:"):
                self.name = str(r[0].split("name:")[1])
                r = f.readline()
                continue
            if r[0].startswith("notes:"):
                self.notes = str(r[0].split("notes:")[1])
                r = f.readline()
                continue
            if r[0].startswith("authors:"):
                self.authors = str(r[0].split("authors:")[1]).split(",")
                r = f.readline()
                continue
            if r[0].startswith("choose_function:"):
                self.choose_function = str(r[0].split("choose_function:")[1])
                r = f.readline()
                continue
            if r[0].startswith("clip_interval:"):
                self.clip_interval = eval(r[0].split("clip_interval:")[1])
                r = f.readline()
                continue
            if r[0].startswith("pair_dist"):
                self.pair_dist = float(r[0].split("pair_dist:")[1])
                r = f.readline()
                continue
            if r[0].startswith("significance_thr"):
                self.significance_thr = float(r[0].split("significance_thr:")[1])
                r = f.readline()
                continue
            if r[0].startswith("iCLIP"):
                fname = r[0].split("iCLIP:")[1]
                self.CLIP.append(fname)
                r = f.readline()
                continue
            if r[0].startswith("CLIP"):
                fname = r[0].split("CLIP:")[1]
                self.CLIP.append(fname)
                r = f.readline()
                continue
            if r[0].startswith("cDNA_thr"):
                self.cDNA_thr = int(r[0].split("cDNA_thr:")[1])
                r = f.readline()
                continue
            if r[0].startswith("presence_thr:"):
                self.presence_thr = float(r[0].split("presence_thr:")[1])
                r = f.readline()
                continue
            if r[0].startswith("control_name:"):
                self.control_name = r[0].split("control_name:")[1]
                r = f.readline()
                continue
            if r[0].startswith("test_name:"):
                self.test_name = r[0].split("test_name:")[1]
                r = f.readline()
                continue
            if r[0].startswith("site_selection"):
                self.site_selection = r[0].split("site_selection:")[1]
                r = f.readline()
                continue
            if r[0].startswith("polya_db:"):
                self.polya_db = r[0].split("polya_db:")[1]
                r = f.readline()
                continue
            if r[0].startswith("poly_type:"):
                self.poly_type = eval(r[0].split("poly_type:")[1])
                r = f.readline()
                continue
            if r[0].startswith("deepbind:"):
                self.deepbind = r[0].split("deepbind:")[1]
                r = f.readline()
                continue
            if r[0].startswith("rnamaps:"):
                self.rnamaps = r[0].split("rnamaps:")[1].split(",")
                r = f.readline()
                continue
            if r[0].startswith("ignore_genes:"):
                self.ignore_genes = r[0].split("ignore_genes:")[1].split(",")
                r = f.readline()
                continue
            if r[0].startswith("exclusive_genes:"):
                self.exclusive_genes = eval(r[0].split("exclusive_genes:")[1])
                r = f.readline()
                continue
            if r[0].startswith("db_type:"):
                self.db_type = r[0].split("db_type:")[1]
                r = f.readline()
                continue
            id = data["id"]
            experiments = data["experiments"].split(",")
            name = data["name"]
            if id.startswith("c"):
                self.control.append((id, experiments, name))
            if id.startswith("t"):
                self.test.append((id, experiments, name))
            r = f.readline()
        f.close()

        # determine comps species
        species = set()
        all_exp = set()
        for (_, exp, _) in self.test+self.control:
            if type(exp)==list:
                for exp_id in exp:
                    all_exp.update(exp)
            else:
                all_exp.update(exp)
        for exp in all_exp:
            lib_id = exp[:exp.rfind("_")]
            exp_id = int(exp.split("_")[-1][1:])
            species.add(apa.annotation.libs[lib_id].experiments[exp_id]["map_to"])
        if len(self.control)>0 and len(self.test)>0:
            assert(len(species)==1)
            self.species = species.pop()

        # finally, load experiment annotation into comps object
        for exp in all_exp:
            lib_id = exp[:exp.rfind("_")]
            exp_id = int(exp.split("_")[-1][1:])
            self.exp_data[exp] = apa.annotation.libs[lib_id].experiments[exp_id]

    def save(self):
        config_file = apa.path.comps_config_filename(self.comps_id)
        f = open(config_file, "wt")
        f.write("\t".join(["id", "experiments", "name"]) + "\n")
        for (id, exp_list, name) in self.control:
            f.write("\t".join(str(x) for x in [id, ",".join(exp_list), name]) + "\n")
        for (id, exp_list, name) in self.test:
            f.write("\t".join(str(x) for x in [id, ",".join(exp_list), name]) + "\n")
        f.write("\n")
        if len(self.CLIP)>0:
            f.write("\n")
            for clip_fname in self.CLIP:
                f.write("CLIP:%s\n" % clip_fname)
            f.write("\n")
        f.write("control_name:%s\n" % self.control_name)
        f.write("test_name:%s\n" % self.test_name)
        f.write("site_selection:%s\n" % self.site_selection)
        f.write("polya_db:%s\n" % self.polya_db)
        f.write("poly_type:%s\n" % json.dumps(self.poly_type))
        f.write("presence_thr:%s\n" % self.presence_thr)
        f.write("cDNA_thr:%s\n" % self.cDNA_thr)
        f.write("\n")
        f.write("authors:%s\n" % ",".join(self.authors))
        f.write("access:%s\n" % ",".join(self.access))
        f.write("owner:%s\n" % ",".join(self.owner))
        f.write("name:%s\n" % self.name)
        f.write("notes:%s\n" % self.notes)
        f.write("\n")
        f.write("analysis_type:%s\n" % self.analysis_type)
        f.write("method:%s\n" % self.method)
        f.write("genome:%s\n" % self.genome)
        f.write("\n")
        f.write("status:%s\n" % self.status)
        f.close()
        return True

def process_comps(comps_id, map_id=1, clean=True):

    if not os.path.exists(os.path.join(apa.path.comps_folder, comps_id, "%s.config" % comps_id)):
        print("%s.config missing, exiting" % (comps_id))
        sys.exit(1)

    # clean
    if clean:
        assert(len(apa.path.comps_folder)>0)
        files = glob.glob(os.path.join(apa.path.comps_folder, comps_id, "*"))
        for f in files:
            if not f.endswith(".config") and not f.startswith("docs") and not f==".htaccess":
                print("removing: %s" % f)
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)

    comps = Comps(comps_id) # study data
    pybio.genomes.load(comps.species)

    flog = open(os.path.join(apa.path.comps_folder, comps_id, "%s.log" % comps_id), "wt")

    # load CLIP data if available
    clip = {}
    for clip_name in comps.CLIP:
        clip[clip_name] = pybio.data.Bedgraph2(os.path.join(apa.path.iCLIP_folder, clip_name))

    if comps.analysis_type=="apa":
        # if there is a polya-db specified in the comparison, load the positions into the filter
        # (strong, weak, less)
        poly_filter = {}
        polydb = pybio.data.Bedgraph()
        for poly_type in comps.poly_type:
            polydb.load(apa.path.polyadb_filename(comps.polya_db, poly_type=poly_type, filetype="bed"), meta=poly_type)
        polydb_annotated = apa.polya.read_polydb(comps.polya_db)

    replicates = []
    expression = {} # keys = c1, c2, c3, t1, t2, t3...items = bedgraph files

    rshort_experiments = {}
    for (comp_id, experiments, comp_name) in comps.control+comps.test:
        key = "%s:%s" % (comp_id, comp_name)
        expression[comp_id] = pybio.data.Bedgraph()
        replicates.append((comp_id, key)) # (short_name, long_name)
        rshort_experiments[comp_id] = experiments
        if comps.analysis_type=="apa":
            for id in experiments:
                print("%s: +%s" % (comp_id, id))
                lib_id = id[:id.rfind("_")]
                exp_id = int(id.split("_")[-1][1:])
                e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=comps.polya_db)
                expression[comp_id].load(e_filename)

    replicates = sorted(replicates, key=lambda item: (item[0][0], int(item[0][1:])))
    print("replicates = ", replicates)

    # save bedgraph files for c1, t1, ...
    beds_folder = os.path.join(apa.path.comps_folder, comps_id, "beds")
    if not os.path.exists(beds_folder):
        os.makedirs(beds_folder)
    dex_folder = os.path.join(apa.path.comps_folder, comps_id, "dex")
    if not os.path.exists(dex_folder):
        os.makedirs(dex_folder)
    if comps.analysis_type=="apa":
        for (rshort, rlong) in replicates:
            fname = os.path.join(beds_folder, "%s.%s.bed" % (comps_id, rshort))
            expression[rshort].save(fname, track_id="%s.%s" % (comps_id, rshort), genome=comps.species)

    # combine all into single file
    # os.system("cat %s/*.bed > %s" % (beds_folder, os.path.join(beds_folder, "%s_all.bed" % comps_id)))

    # global variables
    gsites = {}

    if comps.analysis_type=="apa":

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
            bed_filename = os.path.join(beds_folder, "%s.%s_all.bed.gz" % (comps_id, sample_name))
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
                    (gene_up, gene_id, gene_down, gene_interval, gene_feature) = pybio.genomes.annotate(comps.species, chr, strand, pos)
                    if gene_id==None: # only consider polya sites inside genes
                        continue
                    num_sites_genes += 1
                    num_genes.add(gene_id)
                    sites = gsites.get(gene_id, {})
                    site_data = {"chr":chr, "strand":strand, "pos":pos, "gene_interval":list(gene_interval), "gene_feature":gene_feature} # store position and gene_interval (start, stop, exon/intron), clip binding

                    # get clip data
                    for clip_name in comps.CLIP:
                        # Bedgraph2
                        site_data[clip_name] = clip[clip_name].get_vector("chr"+chr, strand, pos, start=comps.clip_interval[0], stop=comps.clip_interval[1])

                    control_sum = 0
                    test_sum = 0
                    cDNA_sum = 0
                    expression_vector = []
                    for (rshort, _) in replicates:
                        bg = expression[rshort]
                        cDNA = bg.get_value(chr, strand, pos)
                        assert(cDNA>=0)
                        cDNA_sum += cDNA
                        expression_vector.append(cDNA)
                        if rshort.startswith("c"):
                            control_sum += cDNA
                        if rshort.startswith("t"):
                            test_sum += cDNA
                        site_data[rshort] = cDNA

                    # re-added 20170612, not included in the paper version of analysis
                    expression_vector = [1 if cDNA>=comps.cDNA_thr else 0 for cDNA in expression_vector]
                    if sum(expression_vector) < len(expression_vector)/comps.presence_thr:
                        continue

                    # filter lowly expressed sites, paper version
                    if (control_sum<10) and (test_sum<10):
                        continue

                    site_data["cDNA_sum"] = int(cDNA_sum)
                    site_data["site_hex"] = polydb_annotated.get("%s_%s_%s" % (chr, strand, pos), {}).get("PAShex_PASloci_PASindex", "")

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
        matplotlib.use("Agg")
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
        number_moresites = 0
        for e in x:
            y[e] = 0
        for x1 in sorted(num_sites_per_gene.keys()):
            if x1>1:
                number_moresites += num_sites_per_gene[x1]
            if x1<10 and x1!=0:
                y[x1] += num_sites_per_gene[x1]
            if x1>=10:
                y[10] += num_sites_per_gene[x1]
            total_sites += num_sites_per_gene[x1] * x1
            total_genes += num_sites_per_gene[x1]
        y = [y[e] for e in x]
        bar1 = ax.bar(x, y, 0.3, align='center', label='number of polyA sites', color='lightgray')
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        plt.xlim(1-0.25, 10+0.5)
        plt.xlabel('number of polyA sites'); plt.ylabel('number of genes')
        plt.xticks([e for e in x], [1,2,3,4,5,6,7,8,9,">=10"])
        autolabel(bar1)
        plt.title("%s polyA sites annotated to %s genes (%.0f%% genes with APA)" % (format(int(total_sites), ','), format(int(total_genes), ','), number_moresites/float(max(1, total_genes))*100.0))
        plt.tight_layout()
        plt.savefig(os.path.join(apa.path.comps_folder, comps_id, "%s_sites_per_gene.pdf" % comps_id))
        plt.savefig(os.path.join(apa.path.comps_folder, comps_id, "%s_sites_per_gene.png" % comps_id), transparent=True)
        plt.close()
        # =======================

    # expression gene level
    # changed to take counts from htcount of library gene_expression file
    # TODO
    gene_expression_data = {}
    temp_genes = set() # get all genes, do not filter for DGE
    libs = set()
    for (comp_id, experiments, comp_name) in comps.control+comps.test:
        for temp in experiments:
            lib_id = temp[:temp.rfind("_")]
            exp_id = int(temp.split("_")[-1][1:])
            libs.add(lib_id)
    for lib_id in list(libs):
        table_fname = os.path.join(apa.path.data_folder, lib_id, "%s_gene_expression.tab" % lib_id)
        f = open(table_fname, "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            temp_genes.add(r[0])
            for exp_id, exp_value in enumerate(r[2:]):
                gene_expression_data["%s_%s_%s" % (lib_id, exp_id+1, r[0])] = exp_value
            r = f.readline()
        f.close()

    fname = apa.path.comps_expression_filename(comps_id)
    f_genes = open(fname, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype"]
    header += ["sites_position_cDNA", "sites_num", "cDNA_sum"]
    for (rshort, rlong) in replicates:
        header.append(rlong)
    f_genes.write("\t".join(header)+"\n")
    for gene_id in list(temp_genes):
        gene = apa.polya.get_gene(comps.species, gene_id)
        if gene=={}: # gene not found
            continue
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_len = gene_stop-gene_start+1
        assert(gene_len>0)
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)
        row = [chr, strand, gene_locus, gene_id, gene.get("gene_name", ""), gene["gene_biotype"]]
        sites = gsites.get(gene_id, {})
        row.append(sites.keys())
        row.append(len(sites))
        row_comp = []
        cDNA_gtotal = 0
        for (rshort, rlong) in replicates:
            cDNA_rtotal = 0
            for eid in rshort_experiments[rshort]:
                lib_id = eid[:eid.rfind("_")]
                exp_id = int(eid.split("_")[-1][1:])
                ktemp = "%s_%s_%s" % (lib_id, exp_id, gene_id)
                ctemp = int(gene_expression_data.get(ktemp, 0))
                cDNA_rtotal += ctemp
            #cDNA_rtotal = 0
            #for site_pos, es in sites.items():
            #    val = es.get(rshort, 0)
            #    cDNA_rtotal += val
            cDNA_gtotal += cDNA_rtotal
            row_comp.append(cDNA_rtotal)
        row.append(cDNA_gtotal)
        row += row_comp
        f_genes.write("\t".join([str(x) for x in row]) + "\n")
    f_genes.close()

    # expression site level
    gene_sites = {}
    fname = apa.path.comps_expression_filename(comps_id, filetype="sites")
    dex_files = {}
    f_sites = open(fname, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "site_pos", "site_pas", "gene_feature", "gene_interval"]
    for clip_name in comps.CLIP:
        header.append("clip:%s" % clip_name)
    header.append("cDNA_sum")
    for (rshort, rlong) in replicates:
        dex_files[rshort] = open(os.path.join(dex_folder, "dex_%s_%s.tab" % (comps_id, rshort)), "wt")
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
            row_1 = [chr, strand, gene_locus, gene_id, gene["gene_name"], gene["gene_biotype"], site_pos, site_data["site_hex"], site_data["gene_feature"], site_data["gene_interval"]]
            for clip_name in comps.CLIP:
                # TODO
                #row_1.append(site_data[clip_name])
                row_1.append(sum(site_data[clip_name]))
            row_2 = []
            cDNA_sum = 0
            for (rshort, rlong) in replicates:
                val = site_data.get(rshort, 0)
                if len(sites)>1:
                    dex_files[rshort].write("%s_%s%s:%s\t%s\n" % (gene_id, gene["gene_strand"], gene["gene_name"], site_pos, val))
                row_2.append(val)
                cDNA_sum += val
            row = row_1 + [cDNA_sum] + row_2
            f_sites.write("\t".join([str(x) for x in row]) + "\n")
    f_sites.close()
    for (rshort, rlong) in replicates:
        dex_files[rshort].close()

    if clean:

        # differential gene expression with edgeR
        R_file = os.path.join(apa.path.root_folder, "comps", "comps_edgeR.R")
        input_fname = apa.path.comps_expression_filename(comps_id)
        output_fname = os.path.join(apa.path.comps_folder, comps_id, "%s.genes_de.tab" % comps_id)
        command = "R --vanilla --args %s %s %s %s < %s" % (input_fname, output_fname, len(comps.control), len(comps.test), R_file)
        print(command)
        pybio.utils.Cmd(command).run()
        flog.write(command+"\n")
        flog.flush()

        # edgeR normalize gene expression file
        R_file = os.path.join(apa.path.root_folder, "comps", "edgeR_normalize.R")
        input_fname = apa.path.comps_expression_filename(comps_id, filetype="genes")
        output_fname = apa.path.comps_expression_filename(comps_id, filetype="genes_norm")
        command = "R --vanilla --args %s %s < %s" % (input_fname, output_fname, R_file)
        print(command)
        pybio.utils.Cmd(command).run()
        flog.write(command+"\n")
        flog.flush()

        if comps.analysis_type=="apa":
            # replicates cluster analysis
            R_file = os.path.join(apa.path.root_folder, "comps", "comps_cluster.R")
            input_fname = apa.path.comps_expression_filename(comps_id)
            output_fname = os.path.join(apa.path.comps_folder, comps_id, "%s.cluster_genes" % comps_id)
            command = "R --vanilla --args %s %s %s %s %s %s 'analysis: %s, gene expression cluster' < %s" % (input_fname, output_fname, len(comps.control), len(comps.test), comps.cluster_image_w, comps.cluster_image_h, comps_id, R_file)
            print(command)
            pybio.utils.Cmd(command).run()
            flog.write(command+"\n")
            flog.flush()

            # dexseq analysis
            R_file = os.path.join(apa.path.root_folder, "comps", "comps_dex.R")
            output_fname = os.path.join(apa.path.comps_folder, comps_id, "%s" % comps_id)
            command = "R --vanilla --args %s %s %s %s %s < %s" % (dex_folder, output_fname, len(comps.control), len(comps.test), comps_id, R_file)
            print(command)
            pybio.utils.Cmd(command).run()
            flog.write(command+"\n")
            flog.flush()

            # edgeR normalize site expression file
            R_file = os.path.join(apa.path.root_folder, "comps", "edgeR_normalize.R")
            input_fname = apa.path.comps_expression_filename(comps_id, filetype="sites")
            output_fname = apa.path.comps_expression_filename(comps_id, filetype="sites_norm")
            command = "R --vanilla --args %s %s < %s" % (input_fname, output_fname, R_file)
            print(command)
            pybio.utils.Cmd(command).run()
            flog.write(command+"\n")

            # write the pairs_de file
            pairs_de(comps_id, gsites, replicates, polydb)

            # heatmap from pairs_de
            prepare_heatmap_data(comps_id)
            R_file = os.path.join(apa.path.root_folder, "comps", "comps_heatmap_pairsde.R")
            input_fname = apa.path.comps_expression_filename(comps_id)
            command = "R --vanilla --args %s %s %s < %s" % (comps_id, len(comps.control), len(comps.test), R_file)
            print(command)
            pybio.utils.Cmd(command).run()
            flog.write(command+"\n")
            flog.flush()

        flog.close()
        comps.status = "complete"
        comps.save()

def pairs_de(comps_id, gsites, replicates, polydb):
    comps = Comps(comps_id) # study data
    beds_folder = os.path.join(apa.path.comps_folder, comps_id, "beds")
    polydb_annotated = apa.polya.read_polydb(comps.polya_db)
    pairs_filename = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
    f_pairs = open(pairs_filename, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "polyA_sites_in_gene", "proximal_pos", "proximal_pas", "proximal_feature", "proximal_exp", "distal_pos", "distal_pas", "distal_feature", "distal_exp", "s1", "s2"]
    header.append("proximal_control")
    header.append("proximal_control_sum")
    header.append("distal_control")
    header.append("distal_control_sum")
    header.append("proximal_test")
    header.append("proximal_test_sum")
    header.append("distal_test")
    header.append("distal_test_sum")
    header.append("proximal_polyA_type")
    header.append("distal_polyA_type")
    header.append("proximal_fc")
    header.append("distal_fc")
    header.append("proximal_perc_inc")
    header.append("distal_perc_inc")
    header.append("proximal_padj")
    header.append("distal_padj")
    header.append("pair_type")
    header.append("gene_class")

    bg_selected_sites_control = pybio.data.Bedgraph()
    bg_selected_sites_test = pybio.data.Bedgraph()

    # DEX: this completelly invalidates the site selection, and leaves the decision to DEXseq
    if comps.site_selection.upper()=="DEXSEQ":
        selected_sites = dexseq(comps_id, thr=comps.significance_thr)
    else:
        selected_sites = dexseq(comps_id, thr=comps.significance_thr) # for now there is no alternative

    f_pairs.write("\t".join(header)+"\n")

    results = []
    for gene_id, sites in gsites.items():
        gene = apa.polya.get_gene(comps.species, gene_id)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)

        if len(sites)<2:
            continue

        if gene_id not in selected_sites.keys():
            continue

        proximal_pos, distal_pos, proximal_p, distal_p, proximal_fc, distal_fc, gene_class = selected_sites[gene_id]
        proximal_site = sites[proximal_pos]
        distal_site = sites[distal_pos]
        pair_type = apa.polya.annotate_pair(comps.species, chr, strand, proximal_pos, distal_pos)
        s1, s2 = get_s1_s2(gene_id, chr, strand, comps.species, proximal_pos, distal_pos, pair_type)

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
        bg_selected_sites_control.set_value(chr, strand, proximal_pos, sum(proximal_control))
        bg_selected_sites_test.set_value(chr, strand, proximal_pos, sum(proximal_test))
        bg_selected_sites_control.set_value(chr, strand, distal_pos, sum(distal_control))
        bg_selected_sites_test.set_value(chr, strand, distal_pos, sum(distal_test))

        row = [chr, strand, gene_locus, gene_id, gene["gene_name"], gene["gene_biotype"], len(sites)]
        row.append(proximal_pos)
        proximal_hex = polydb_annotated.get("%s_%s_%s" % (chr, strand, proximal_pos), {}).get("PAShex_PASloci_PASindex", "")
        row.append(proximal_hex)
        row.append(proximal_site["gene_feature"])
        row.append(sum(proximal_test+proximal_control))
        row.append(distal_pos)
        distal_hex = polydb_annotated.get("%s_%s_%s" % (chr, strand, distal_pos), {}).get("PAShex_PASloci_PASindex", "")
        row.append(distal_hex)
        row.append(distal_site["gene_feature"])
        row.append(sum(distal_test+distal_control))

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

        # strong, weak, noclass...
        #proximal_class = polydb.get_value(chr, strand, proximal_pos, db="meta")
        #distal_class = polydb.get_value(chr, strand, distal_pos, db="meta")
        proximal_class = polydb_annotated.get("%s_%s_%s" % (chr, strand, proximal_pos), {}).get("pas_type", "")
        distal_class = polydb_annotated.get("%s_%s_%s" % (chr, strand, distal_pos), {}).get("pas_type", "")
        row.append(proximal_class)
        row.append(distal_class)

        row.append("%.5f" % proximal_fc)
        row.append("%.5f" % distal_fc)

        row.append("%.2f" % (sum(proximal_control)/float(sum(proximal_control)+sum(proximal_test))))
        row.append("%.2f" % (sum(distal_control)/float(sum(distal_control)+sum(distal_test))))
        row.append("%.5f" % proximal_p)
        row.append("%.5f" % distal_p)

        pair_type = apa.polya.annotate_pair(comps.species, chr, strand, proximal_site["pos"], distal_site["pos"])
        row.append(pair_type)
        row.append(gene_class)
        results.append(row)

    results = sorted(results, key=lambda x: x[-1], reverse=True)
    results = sorted(results, key=lambda x: x[-2], reverse=True)

    # write in this order
    for rt in ["repressed", "enhanced", "control_up", "control_down"]:
        for pt in ["same", "skipped", "composite", "other"]:
            for row in results:
                if row[-2]==pt and row[-1]==rt:
                    f_pairs.write("\t".join([str(x) for x in row]) + "\n")
#    for row in results:
#        f_pairs.write("\Z" + "\t".join([str(x) for x in row]) + "\n")
    f_pairs.close()

    # save selected sites bedGraphs
    bg_selected_sites_control_fname = os.path.join(beds_folder, "%s_control_selected.bed.gz" % comps_id)
    bg_selected_sites_test_fname = os.path.join(beds_folder, "%s_test_selected.bed.gz" % comps_id)
    bg_selected_sites_control.save(bg_selected_sites_control_fname, track_id="%s_control_selected" % comps_id)
    bg_selected_sites_test.save(bg_selected_sites_test_fname, track_id="%s_test_selected" % comps_id)

"""
Get splice sites related to polyA sites.
https://docs.google.com/drawings/d/1_m4iZ1c9YwKI-NOWMSCSdg2IzEGedj-NaMTIlHqirc0
"""
def get_s1_s2(gene_id, chr, strand, genome, proximal_pos, distal_pos, pair_type):
    pybio.genomes.load(genome)
    gene = apa.polya.get_gene(genome, gene_id)
    _, _, _, gene_interval, gene_feature = pybio.genomes.annotate(genome, chr, strand, proximal_pos)
    proximal_site = {"pos":proximal_pos, "gene_interval":list(gene_interval), "gene_feature":gene_feature}
    _, _, _, gene_interval, gene_feature = pybio.genomes.annotate(genome, chr, strand, distal_pos)
    distal_site = {"pos":distal_pos, "gene_interval":list(gene_interval), "gene_feature":gene_feature}

    s1, s2 = None, None

    if pair_type=="same":
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

"""
Draw hisogram of proximal and distal polyA site distances per gene.
"""
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

"""
Use DEXSeq to select 2 polyA sites per gene and classify the gene into repressed / enhanced / control_up / control_down category
"""
def dexseq(comps_id, thr=0.05):
    comps = Comps(comps_id) # study data
    f = open(os.path.join(apa.path.comps_folder, comps_id, "%s.dex.tab" % comps_id), "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    results = {}
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        gene_id = data["groupID"]
        L = results.get(gene_id, [])
        gene_exp = 0
        for el in r[11:]:
            gene_exp += int(el)
        strand = data["groupID"].split("_")[1][0]
        try:
            pos = int(data["featureID"][1:])
        except:
            print("Could not recover position from: %s" % data["featureID"])
            r = f.readline()
            continue
        if data["padj"]=="NA" or data["log2fold_test_control"]=="NA":
            r = f.readline()
            continue
        # our DEXSeq script reports fold-changes with log2_test_control orientation, we reverse the sign to report in log2_control_test orientation
        row = {"padj": float(data["padj"]), "fc": -float(data["log2fold_test_control"]), "featureID": data["featureID"], "groupID": data["groupID"], "gene_exp": gene_exp, "strand": strand, "pos": pos}
        L.append(row)
        results[gene_id] = L
        r = f.readline()
    f.close()

    dex_results = {}
    repressed = 0; enhanced = 0; c_up = 0; c_down = 0
    for gene_id, L in results.items():
        gene_class = None
        pair_type = None
        gid = gene_id.split("_")[0]
        # assert(len(L)>1) # no longer valid, since if the expression accross replicates is too low, DEXseq returns NA for certain sites and it could be that certain genes have only 1 polyA site
        if len(L)<2:
            continue
        sig_sites = [x for x in L if x["padj"]<=thr]
        control_sites = [x for x in L if x["padj"]>thr]

        if len(sig_sites)>1:
            pairs = list(itertools.combinations(sig_sites, 2))
            all_pairs = []
            for L1, L2 in pairs:
                if L1["fc"]*L2["fc"]<0 and abs(L1["pos"]-L2["pos"])>comps.pair_dist: # direction opposite and pair distant enough? consider pair
                    all_pairs.append((abs(L1["fc"]-L2["fc"]), L1, L2))
            #all_pairs.sort(reverse=True)
            all_pairs.sort(key=lambda x: x[0], reverse=True) # python 3
            if len(all_pairs)>0:
                site1, site2 = all_pairs[0][1], all_pairs[0][2]
                pair_type = "reg"
        elif len(control_sites)>1:
            control_sites.sort(key=lambda x: x["gene_exp"], reverse=True)
            all_pairs = []
            for L1, L2 in zip(control_sites, control_sites[1:]): # 1,2; 2,3; 3,4; 4,5;...
                if abs(L1["pos"]-L2["pos"])>comps.pair_dist:
                    all_pairs.append((L1["gene_exp"]+L2["gene_exp"], L1, L2))
            #all_pairs.sort(reverse=True) # just to check, because they are already sorted by highest expression
            all_pairs.sort(key=lambda x: x[0], reverse=True) # python 3
            if len(all_pairs)>0:
                site1, site2 = all_pairs[0][1], all_pairs[0][2]
                pair_type = "control"
        elif len(sig_sites)==1:
            all_pairs = []
            site1 = sig_sites[0]
            control_sites.sort(key=lambda x: x["gene_exp"], reverse=True)
            for potential_control in control_sites:
                if abs(site1["pos"]-potential_control["pos"])>comps.pair_dist:
                    all_pairs.append((site1, potential_control))
            if len(all_pairs)>0:
                site1, site2 = all_pairs[0][0], all_pairs[0][1]
                pair_type = "control"

        # determine proximal and distal
        if pair_type in ["reg", "control"]:
            if site1["pos"]<site2["pos"]:
                proximal = site1
                distal = site2
            else:
                proximal = site2
                distal = site1
            strand = site1["strand"]

            if strand=="-":
                proximal, distal = distal, proximal

        # gene regulation
        if pair_type=="control":
            if proximal["fc"]<0:
                c_up += 1
                gene_class = "control_up"
            else:
                c_down += 1
                gene_class = "control_down"
        elif pair_type=="reg":
            if proximal["fc"]>0:
                gene_class = "enhanced"
                enhanced += 1
            else:
                gene_class = "repressed"
                repressed += 1

        if gene_class!=None:
            proximal_pos = proximal["pos"]
            distal_pos = distal["pos"]
            proximal_p = proximal["padj"]
            distal_p = distal["padj"]
            distance = abs(proximal_pos-distal_pos)
            proximal_fc = proximal["fc"]
            distal_fc = distal["fc"]
            dex_results[gid] = (proximal_pos, distal_pos, proximal_p, distal_p, proximal_fc, distal_fc, gene_class)
    return dex_results

"""
Plotly of fold changes from pairs_de table
"""
def apa_plot(comps_id):
    pairs_filename = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
    html_filename = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.html" % comps_id)
    plot_data = {"enhanced" : {"x":[], "y":[], "gene_id":[]}, "repressed" : {"x":[], "y":[], "gene_id":[]}, "control_up" : {"x":[], "y":[], "gene_id":[]}, "control_down" : {"x":[], "y":[], "gene_id":[]}}
    f = open(pairs_filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        plot_data[data["gene_class"]]["x"].append(float(data["proximal_fc"]))
        plot_data[data["gene_class"]]["y"].append(float(data["distal_fc"]))
        plot_data[data["gene_class"]]["gene_id"].append("%s: %s" % (data["gene_id"], data["gene_name"]))
        r = f.readline()
    f.close()

    html_text = """
<html>
<head>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>

<style>
  div.plotly-notifier {{
    visibility: hidden;
  }}
</style>

<body>

  <div id="myDiv" style="width: 100%%; 100%%;"></div>
  <script>

      var data_e = {{
        x : {enhanced_x},
        y : {enhanced_y},
        mode: 'markers',
        type: 'scatter',
        name: 'Enhanced {enhanced_num}',
        text: {enhanced_gene_id},
        marker: {{ size: 8, color: "#ff6961" }}
      }};

      var data_r = {{
        x : {repressed_x},
        y : {repressed_y},
        mode: 'markers',
        type: 'scatter',
        name: 'Repressed {repressed_num}',
        text: {repressed_gene_id},
        marker: {{ size: 8, color: "#779ecb" }}
      }};

      var data_cup = {{
        x : {control_up_x},
        y : {control_up_y},
        mode: 'markers',
        type: 'scatter',
        name: 'Controls up {control_up_num}',
        text: {control_up_gene_id},
        marker: {{ size: 8, color: "#f1f1f1" }}
      }};

      var data_cdown = {{
        x : {control_down_x},
        y : {control_down_y},
        mode: 'markers',
        type: 'scatter',
        name: 'Controls up {control_down_num}',
        text: {control_down_gene_id},
        marker: {{ size: 8, color: "#f1f1f1" }}
      }};

      var data = [ data_cup, data_cdown, data_e, data_r];

        var layout = {{
          margin: {{
            l: 50,
            r: 50,
            b: 50,
            t: 20,
            pad: 4
          }},
          xaxis: {{
            title: "proximal polyA site (fold-change)"
          }},
          yaxis: {{
            title: "distal polyA site (fold-change)",
          }},
          title: '{title}',
          hovermode:'closest',
          font: {{
            family: 'Arial',
            size: 10,
            color: '#7f7f7f'
          }}
      }};

        Plotly.newPlot('myDiv', data, layout, {{displayModeBar: false}});
    </script>
    </body>
"""

    final_text = html_text.format(title=comps_id, enhanced_x=plot_data["enhanced"]["x"], enhanced_y=plot_data["enhanced"]["y"], enhanced_gene_id=plot_data["enhanced"]["gene_id"], enhanced_num=len(plot_data["enhanced"]["gene_id"]), repressed_x=plot_data["repressed"]["x"], repressed_y=plot_data["repressed"]["y"], repressed_gene_id=plot_data["repressed"]["gene_id"], repressed_num=len(plot_data["repressed"]["gene_id"]), control_up_x=plot_data["control_up"]["x"], control_up_y=plot_data["control_up"]["y"], control_up_gene_id=plot_data["control_up"]["gene_id"], control_up_num=len(plot_data["control_up"]["gene_id"]), control_down_x=plot_data["control_down"]["x"], control_down_y=plot_data["control_down"]["y"], control_down_gene_id=plot_data["control_down"]["gene_id"], control_down_num=len(plot_data["control_down"]["gene_id"]))
    open(html_filename, "wt").write(final_text)

    # also store a PDF version of the plot
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
    matplotlib.rcParams['axes.labelsize'] = 8
    matplotlib.rcParams['axes.titlesize'] = 8
    matplotlib.rcParams['xtick.labelsize'] = 8
    matplotlib.rcParams['ytick.labelsize'] = 8
    matplotlib.rcParams['legend.fontsize'] = 8
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.4
    matplotlib.rcParams['legend.frameon'] = 'False'

    for c_pair_type in ["same", "composite", "skipped", "combined"]:
        plot_data = {"enhanced" : {"x":[], "y":[], "gene_id":[]}, "repressed" : {"x":[], "y":[], "gene_id":[]}, "control_up" : {"x":[], "y":[], "gene_id":[]}, "control_down" : {"x":[], "y":[], "gene_id":[]}}
        f = open(pairs_filename, "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            if data["pair_type"]==c_pair_type or c_pair_type=="combined":
                plot_data[data["gene_class"]]["x"].append(float(data["proximal_fc"]))
                plot_data[data["gene_class"]]["y"].append(float(data["distal_fc"]))
                plot_data[data["gene_class"]]["gene_id"].append("%s: %s" % (data["gene_id"], data["gene_name"]))
            r = f.readline()
        f.close()
        fig, ax = plt.subplots(1, 1, figsize=(12, 3))
        marksize = 20
        markalpha = 0.9
        ax.scatter(plot_data["control_up"]["x"], plot_data["control_up"]["y"], color='#f1f1f1', s=marksize, edgecolor='none', alpha=markalpha)
        ax.scatter(plot_data["control_down"]["x"], plot_data["control_down"]["y"], color='#f1f1f1', s=marksize, edgecolor='none', alpha=markalpha)
        ax.scatter(plot_data["enhanced"]["x"], plot_data["enhanced"]["y"], color='#ff6961', s=marksize, edgecolor='none', alpha=markalpha)
        ax.scatter(plot_data["repressed"]["x"], plot_data["repressed"]["y"], color='#779ecb', s=marksize, edgecolor='none', alpha=markalpha)
        plt.title("%s '%s' APA switches %s repressed, %s enhanced, %s control_up, %s control_down" % (c_pair_type, comps_id, len(plot_data["repressed"]["x"]), len(plot_data["enhanced"]["x"]), len(plot_data["control_up"]["x"]), len(plot_data["control_down"]["x"])))
        plt.tight_layout()
        plt.savefig(os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.%s.pdf" % (comps_id, c_pair_type)))
        plt.close()

def prepare_heatmap_data(comps_id):
    f = open(os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id), "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    ntest = None
    ncontrol = None
    ngenes = 0
    results = [] # already computed PSI, stored to heatmap.tab
    results_complete = [] # all numbers, stored to heatmap_complete.tab

    comps = apa.comps.Comps(comps_id)

    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        data = dict(zip(header, r))
        pcontrol = data["proximal_control"].split(";")
        pcontrol = [float(x) for x in pcontrol]
        ptest = data["proximal_test"].split(";")
        ptest = [float(x) for x in ptest]
        dcontrol = data["distal_control"].split(";")
        dcontrol = [float(x) for x in dcontrol]
        dtest = data["distal_test"].split(";")
        dtest = [float(x) for x in dtest]

        assert(len(dtest)==len(ptest))
        assert(len(dcontrol)==len(pcontrol))

        if ncontrol==None:
            ncontrol = len(pcontrol)
        if ntest==None:
            ntest = len(ptest)

        row = ["%s %s" % (data["gene_name"], data["gene_class"])]
        row_complete = ["%s %s" % (data["gene_name"], data["gene_class"])]
        for i in range(0, ncontrol):
            dv = pcontrol[i]+dcontrol[i]
            if dv==0:
                dv = 1
            row.append(int(pcontrol[i]/dv*100))
            row_complete.append("%s;%s" % (int(pcontrol[i]), int(dcontrol[i])))
        for i in range(0, ntest):
            dv = ptest[i]+dtest[i]
            if dv==0:
                dv = 1
            row.append(int(ptest[i]/dv*100))
            row_complete.append("%s;%s" % (int(ptest[i]), int(dtest[i])))
        average_control = sum(row[1:1+ncontrol])/ncontrol
        average_test = sum(row[ncontrol:])/ntest
        dPSI = abs(average_control - average_test)

        if data["gene_class"] in ["repressed", "enhanced"]: # and data["pair_type"] in ["tandem"]:
            results.append((dPSI, row))
            results_complete.append((dPSI, row_complete))
            ngenes += 1

        r = f.readline()
    f.close()

    fout = open(os.path.join(apa.path.comps_folder, comps_id, "%s.heatmap.tab" % comps_id), "wt")
    header2 = ["gene_name"]
    for i in range(0, ncontrol):
        header2.append(comps.control[i][2])
    for i in range(0, ntest):
        header2.append(comps.test[i][2])
    fout.write("\t".join(header2)+"\n")
    results.sort(reverse=True)
    for (dPSI, row) in results:
        fout.write("\t".join(str(x) for x in row) + "\n")
    fout.close()

    fout_complete = open(os.path.join(apa.path.comps_folder, comps_id, "%s.complete_heatmap.tab" % comps_id), "wt")
    header2 = ["gene_name"]
    for i in range(0, ncontrol):
        header2.append(comps.control[i][2])
    for i in range(0, ntest):
        header2.append(comps.test[i][2])
    fout_complete.write("\t".join(header2)+"\n")
    results.sort(reverse=True)
    for (dPSI, row) in results_complete:
        fout_complete.write("\t".join(str(x) for x in row) + "\n")
    fout_complete.close()
