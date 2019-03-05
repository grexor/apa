import pysam
import sys
import apa
import os
import pybio
import time
import glob
import gzip
import numpy as np
import shutil
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['axes.titlesize'] = 10
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['legend.fontsize'] = 9

def get_species(poly_id):
    if poly_id in ["hg19_tian", "hg19_derti"]:
        return "hg19"
    if poly_id in pybio.genomes.genomes_list():
        species = poly_id
    else: # ? poly_id is a pool of experiments? read config file
        species = set()
        f = open(os.path.join(apa.path.polya_folder, "%s.config" % poly_id), "rt")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "")
            lib_id = "_".join(r.split("_")[:-1])
            exp_id = int(r.split("_")[-1][1:])
            map_to = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
            species.add(map_to)
            r = f.readline()
        f.close()
        assert(len(species)==1)
        species = species.pop()
    return species

def make_config(lib_id):
    f = open(os.path.join(apa.path.polya_folder, "%s.config" % lib_id), "wt")
    lib = apa.annotation.libs[lib_id]
    for exp_id, exp_data in lib.experiments.items():
        f.write("%s_e%s\n" % (lib_id, exp_id))
    f.close()

def process(poly_id, map_id=1, min_distance=25):
    """
    Creates polyA database.
    It loads multiple experiments (T files). Scales each experiment individually and adds up the raw counts and the scaled values.
    Also, each experiment is identified with an id (filename). All this is done by pybio.data.Bedgraph.load()
    cpm: counts per milion: when scaled, each experiment position if tested: scale vaue >= cpm, and file id is added to the "support" database accordingly
    After all T files are combined and read into the master bed file, clustering is performed on raw data, scale data, support data and meta data.

    bed.save(polyadb_filename+".%02gmincDNA.%03gcpm" % (min_cDNA, cpm), db_filter="support", db_save="raw", min_cDNA=min_cDNA)

    This saves all positions where the minimal number of experiments that contributed to the position >= min_cDNA. This is the filtering, but
    raw data (sum of all experiments) is saved.
    """
    species = get_species(poly_id)

    experiments = []
    if poly_id in pybio.genomes.genomes_list():
        # is the poly_id a genome assembly name?
        for lib_id, lib_data in apa.annotation.libs.items():
            for exp_id, exp_data in lib_data.experiments.items():
                if poly_id!=exp_data["map_to"]:
                    continue
                experiments.append((lib_id, exp_id))
    else:
        # is the poly_id a list of experiments (=pool) filename?
        f = open(os.path.join(apa.path.polya_folder, "%s.config" % poly_id), "rt")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "")
            lib_id = "_".join(r.split("_")[:-1])
            exp_id = int(r.split("_")[-1][1:])
            experiments.append((lib_id, exp_id))
            map_to = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
            r = f.readline()
        f.close()

    bed = pybio.data.Bedgraph()
    num_read = 0
    for lib_id, exp_id in experiments:
        num_read += 1
        meta = ["%s_e%s" % (lib_id, exp_id)] # list of metadata
        t_filename = apa.path.t_filename(lib_id, exp_id, map_id=map_id)
        print t_filename
        if os.path.exists(t_filename):
            bed.load(t_filename, meta=meta)
            print "%s: %s %s %s %s %.2fM" % (num_read, lib_id, exp_id, poly_id, os.path.exists(t_filename), bed.total_raw/1000000.0)

    #bed.filter(min_distance=25) # Gregor: alternative: -25..25 (201702 test)
    bed.filter(min_distance=min_distance)
    bed.save(apa.path.polyadb_filename(poly_id, filetype="temp"), db_save="raw")

    annotate(poly_id)

def get_gene(species, gid):
    # only used by apa.sition
    # return only one number of gene_start and gene_stop
    # if no gene, return {}
    pybio.genomes.load(species) # make sure the annotation is loaded
    gene = pybio.genomes.genes.get(species, {}).get(gid, {})
    if gene!={}:
        if type(gene["gene_start"])==list:
            gene["gene_start"] = gene["gene_start"][0]
        if type(gene["gene_stop"])==list:
            gene["gene_stop"] = gene["gene_stop"][-1]
    return gene

def annotate(poly_id):
    species = get_species(poly_id)
    interval_types = {"3":"3utr", "5": "5utr", "o":"orf", "i":"intron"}
    polyadb_temp = apa.path.polyadb_filename(poly_id, filetype="temp")
    polyadb_tab = apa.path.polyadb_filename(poly_id, filetype="tab")
    polyadb_bed = apa.path.polyadb_filename(poly_id, filetype="bed")
    polyadb_fasta = apa.path.polyadb_filename(poly_id, filetype="fasta")

    db = {}
    db_values = []
    gene_values = {}
    print polyadb_temp
    f = gzip.open(polyadb_temp)
    r = f.readline()
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        chr = r[0].replace("chr", "")
        pos = int(r[1])
        cDNA = float(r[-1])
        strand = "+" if cDNA>=0 else "-"
        gid_up, gid, gid_down, gid_interval, _ = pybio.genomes.annotate(species, chr, strand, pos)
        key = "%s:%s:%s" % (chr, strand, pos)
        db[key] = (cDNA, gid, gid_interval)
        db_values.append(abs(cDNA))
        if gid!=None:
            gene_values.setdefault(gid, []).append(abs(cDNA))
        r = f.readline()
    f.close()

    # at least cDNA 10 to accept polyA database position
    cDNA_filter = {"hg19_derti":0, "hg19_tian":0}.get(poly_id, 10)
    accepted = {}
    for key, (cDNA, gid, gid_interval) in db.items():
        if abs(cDNA)>cDNA_filter:
            accepted[key] = 1

    ftab = gzip.open(polyadb_tab, "wb")
    ftab.write("\t".join(["chr", "strand", "pos", "gene_id", "gene_name", "interval", "cDNA", "pas_type", "pas_loci", "cs_loci", "seq_-100_100"]) + "\n")
    fbed = gzip.open(polyadb_bed, "wb")
    fbed.write("track type=bedGraph name=\"%s\" description=\"%s\" altColor=\"200,120,59\" color=\"120,101,172\" maxHeightPixels=\"100:50:0\" visibility=\"full\" priority=\"20\"" % (poly_id, poly_id))

    # store upstream sequences for polyar to classify weak/strong/other poly-A sites
    ffasta = open(polyadb_fasta, "wt")
    f = gzip.open(polyadb_temp)
    r = f.readline()
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        chr = r[0].replace("chr", "")
        pos = int(r[1])
        cDNA = float(r[-1])
        strand = "+" if cDNA>=0 else "-"
        # get upstream sequence
        seq = pybio.genomes.seq(species, chr, strand, pos, start=-100, stop=100)
        gid_up, gid, gid_down, gid_interval, _ = pybio.genomes.annotate(species, chr, strand, pos)
        key = "%s:%s:%s" % (chr, strand, pos)
        if accepted.get(key, None)!=None:
            if gid==None:
                row = [chr, strand, pos, "", "", "", cDNA, "", "", "", seq]
            else:
                gene = get_gene(species, gid)
                interval = "%s:%s:%s" % (str(gid_interval[0]), str(gid_interval[1]), interval_types[gid_interval[2]])
                row = [chr, strand, pos, gid, gene["gene_name"], interval, cDNA, "", "", "", seq]
            ftab.write("\t".join(str(e) for e in row)+"\n")
            fbed.write("\t".join(r)+"\n")
            ffasta.write(">%s\n%s\n" % ("%s%s:%s" % (strand, chr, pos), seq))
        r = f.readline()
    f.close()
    ftab.close()
    fbed.close()
    ffasta.close()
    os.remove(polyadb_temp)
    classify_polya(poly_id)
    polyadb_class_histogram(poly_id)

def classify_polya(poly_id):
    """
    Create a fasta file from the -100, 100 sequence around detected poly-A sites.
    Use polyar to classify those sites.
    Update the polya_db.tab file to include site type (strong, weak, less, "").
    Empty means no CS detected;
    Also add the location of the CS and PAS signal.
    """

    polyadb_tab = apa.path.polyadb_filename(poly_id, filetype="tab")
    polyadb_fasta = apa.path.polyadb_filename(poly_id, filetype="fasta")
    shutil.copy2(polyadb_fasta, os.path.join(os.getenv("HOME"), "software/polyar/"))
    os.chdir(os.path.join(os.getenv("HOME"), "software/polyar/"))

    # run polyar
    os.system("java polyar -i %s -o %s.res" % (polyadb_fasta, poly_id))

    polyar_results = {}
    f = open("%s.res" % poly_id, "rt")
    r = f.readline()
    text = ""
    while r:
    	text = text + r
    	r = f.readline()
    f.close()
    text = text.split("Query Name")
    for line in text[1:]:
            line = line.split("\n")
            id = line[0][4:]
            strand = id[0]
            chr = id.split(":")[0][1:]
            pos = int(id.split(":")[1])
            pas_type = ""
            cs = ""
            pas = ""
            if line[2]!="":
                    if line[2].find("strong")!=-1:
                            pas_type = "strong"
                    elif line[2].find("weak")!=-1:
                            pas_type = "weak"
                    elif line[2].find("less")!=-1:
                            pas_type = "less"
                    cs = int(line[3].replace("\t", "").split("CS:")[1][:5]) - 100 # because sequence is -100, 100, 100 = 0
                    if line[3].replace("\t", "").find("PAS")!=-1:
                        pas = int(line[3].replace("\t", "").split("PAS:")[1][:5]) - 100 # because sequence is -100, 100, 100 = 0
            polyar_results[(chr, strand, pos)] = (pas_type, pas, cs)

    f = gzip.open(polyadb_tab)
    files = {}
    f_tab = gzip.open("%s.tab.gz" % poly_id, "wb")
    for poly_type in ["tab", "strong", "weak", "less", "noclass", "p31"]:
        files[poly_type] = gzip.open(apa.path.polyadb_filename(poly_id, poly_type=poly_type, filetype="bed"), "wb")
    header = f.readline()
    f_tab.write(header)
    r = f.readline()
    while r:
    	r = r.replace("\n", "").replace("\r", "").split("\t")
    	chr, strand, pos, gene_id, gene_name, interval, cDNA, seq = r[0], r[1], int(r[2]), r[3], r[4], r[5], float(r[6]), r[10]
        pas_type, pas_offset, cs_offset = polyar_results[(chr, strand, pos)]
        if pas_type=="":
            pas_type="noclass"
        r[-4], r[-3], r[-2] = pas_type, pas_offset, cs_offset
    	f_tab.write("\t".join([str(e) for e in r])+"\n")
        bed_row = [chr, pos, pos+1, cDNA]
        files[pas_type].write("\t".join(str(e) for e in bed_row) + "\n")
    	r = f.readline()
    f.close()
    f_tab.close()
    for f in files.values():
        f.close()

    os.system("cp %s.tab.gz %s" % (poly_id, polyadb_tab))
    #os.system("rm %s.res" % poly_id)
    #os.system("rm %s.fasta" % poly_id)
    #os.system("rm %s.tab" % poly_id)
    return

def polyadb_class_histogram(poly_id):

    polyadb_tab = apa.path.polyadb_filename(poly_id, filetype="tab")
    f = gzip.open(polyadb_tab)
    y = {"strong":[], "weak":[], "less":[], "noclass":[]}
    header = f.readline().replace("\n", "").replace("\r", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        data = dict(zip(header, r))
        if data["pas_type"] not in ["noclass"]:
            y[data["pas_type"]].append(int(data["cs_loci"]))
        else:
            y[data["pas_type"]].append(0) # just for counting by len, no drawing
        r = f.readline()
    f.close()

    import matplotlib
    import numpy as np
    import pylab as P
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tkr
    # styling
    matplotlib.rcParams['axes.labelsize'] = 15
    matplotlib.rcParams['axes.titlesize'] = 15
    matplotlib.rcParams['xtick.labelsize'] = 13
    matplotlib.rcParams['ytick.labelsize'] = 13
    matplotlib.rcParams['legend.fontsize'] = 13
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.1
    matplotlib.rcParams['legend.frameon'] = 'False'
    from matplotlib import ticker

    def func(x, pos):  # formatter function takes tick label and tick position
       s = '{:0,d}'.format(int(x))
       return s

    P.figure(figsize=(12,3))
    n, bins, patches = P.hist([y["strong"], y["weak"], y["less"], []], bins=np.arange(-30, 30+1)-0.5, color=["#FFAEAE", 'lightgreen', "lightblue", "#f1f1f1"], label=["strong (%s)" % format(int(len(y["strong"])), ","), "weak (%s)" % format(int(len(y["weak"])), ","), "PAS-less (%s)" % format(int(len(y["less"])), ","), "PAS-less (no predicted CS) (%s)" % format(int(len(y["noclass"])), ",")], histtype="barstacked", edgecolor="none")
    P.title("distribution of predicted polyAR cleavage sites around %s cleavage sites" % poly_id)
    P.ylim(bottom=0)
    P.xlabel("distance from cleavage site [nt]")
    P.ylabel("number of sites")
    P.legend()

    a_format = tkr.FuncFormatter(func)
    axes = plt.gca()
    axes.xaxis.set_major_formatter(a_format)
    axes.yaxis.set_major_formatter(a_format)
    axes.spines['bottom'].set_alpha(0.5)
    axes.spines['top'].set_alpha(0.5)
    axes.spines['right'].set_alpha(0.5)
    axes.spines['left'].set_alpha(0.5)
    axes.set_xlim(-20, 20)
    polyadb_image = apa.path.polyadb_filename(poly_id, filetype="polyar_pdf")
    P.tight_layout()
    P.savefig(polyadb_image)

def annotate_pair(species, chr, strand, pos1, pos2, extension=5000):
    # keep positions in positive orientation
    if pos1>pos2:
        pos1, pos2 = pos2, pos1
    _, gid1, _, interval1, _ = pybio.genomes.annotate(species, chr, strand, pos1, extension=extension)
    _, gid2, _, interval2, _ = pybio.genomes.annotate(species, chr, strand, pos2, extension=extension)

    assert(gid1==gid2) # pair needs to be in the same gene
    # intervals are always in the positive orientation: index0, index1, index2...regardless of strand
    # example: - strand: [[100,200,'o']=3'UTR, [200,300,'i'], [300,400,'o']=5'UTR]
    gene_intervals = get_gene(species, gid1)["gene_intervals"]
    index1 = gene_intervals.index(list(interval1))
    index2 = gene_intervals.index(list(interval2))
    upint = None
    if strand=="+" and index1>0:
        upint = gene_intervals[index1-1][2]
    if strand=="-" and index2<len(gene_intervals)-1:
        upint = gene_intervals[index2+1][2]
    intervals = [gene_intervals[i] for i in range(index1, index2+1)]
    itypes = [i[2] for i in intervals]
    itypes = "".join(itypes) # e.g.: "5oioioi3", in this case site1 would be in 5, site2 in 3 (at extremes of itypes)
    if strand=="-":
        itypes = itypes[::-1] # reverse string
    #if itypes.find("io")==-1 and itypes.find("oi")==-1:
    if itypes=="o": # tandem sites need to be in the same exon
        return "same"
    if itypes.count("io")!=-1:
        # check upstream interval of upstream site
        if upint==None:
            return "skipped"
        if "%s%s" % (upint, itypes[0])=="oi":
            return "composite"
        return "skipped"
    return "other"

def pas_db(poly_id, map_id=1):
    """
    Creates PAS database.
    """
    species = get_species(poly_id)
    min_distance = 60

    experiments = []
    if poly_id in pybio.genomes.genomes_list():
        # is the poly_id a genome assembly name?
        for lib_id, lib_data in apa.annotation.libs.items():
            for exp_id, exp_data in lib_data.experiments.items():
                if poly_id!=exp_data["map_to"]:
                    continue
                experiments.append((lib_id, exp_id))
    else:
        # is the poly_id a list of experiments (=pool) filename?
        f = open(os.path.join(apa.path.polya_folder, "%s.config" % poly_id), "rt")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "")
            lib_id = "_".join(r.split("_")[:2])
            exp_id = int(r.split("_")[-1][1:])
            experiments.append((lib_id, exp_id))
            map_to = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
            r = f.readline()
        f.close()

    # debug reasons
    #experiments = experiments[:2]

    b = pybio.data.Bedgraph()
    for (lib_id, exp_id) in experiments:
        fname = apa.path.t_filename(lib_id, exp_id, map_id=map_id)
        b.load(fname)
    name_t = os.path.join(apa.path.polya_folder, "%s.data_t" % poly_id)
    b.save(name_t+".bed.gz")

    b = pybio.data.Bedgraph()
    for (lib_id, exp_id) in experiments:
        fname = apa.path.r_filename(lib_id, exp_id, map_id=map_id)
        b.load(fname)
    name_r = os.path.join(apa.path.polya_folder, "%s.data_r" % poly_id)
    b.save(name_r+".bed.gz")

    def write_seqs(tab_filename):
        pybio.genomes.load("hg19")
        f = gzip.open(tab_filename, "wb")
        f.write("\t".join(["chr", "strand", "pos", "cDNA", "s_-100_100"]) + "\n")
        for chr, chr_data in b.raw.items():
            for strand, pos_data in chr_data.items():
                for pos, val in pos_data.items():
                    s = pybio.genomes.seq("hg19", chr, strand, pos, start=-100, stop=100)
                    f.write("\t".join(str(el) for el in [chr, strand, pos, val, s])+"\n")
        f.close()

    pybio.genomes.load(species)
    write_seqs(name_t+".tab.gz")
    write_seqs(name_r+".tab.gz")

    db = pybio.data.Bedgraph()
    db_pas = {}

    pas_signals = ["AATAAA", "ATTAAA", "AAATAA", "CAATAA", "AGTAAA", "TTAATA", "TATAAA", "TGAATA", "CATAAA", "GTAATA", "CTAATA", "GATAAA"]

    f = gzip.open(name_t+".tab.gz")
    header = f.readline()
    r = f.readline()
    all_seqs = 0
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        all_seqs += 1
        seq = r[-1][63:87] # upstream 100 nt
        cDNA = int(r[-2])
        pos = int(r[-3])
        strand = r[-4]
        chr = r[-5]
        for pas in pas_signals:
            index = seq.find(pas)
            if index!=-1:
                start = -(100-(63+index))
                pas_check = pybio.genomes.seq("hg19", chr, strand, pos, start=start, stop=start+5)
                assert(pas==pas_check)
                db.set_value(chr, strand, pos, db.get_value(chr, strand, pos)+cDNA) # increase the value by cDNA at this position
                db_pas.setdefault((chr, strand, pos), set()).add(pas)
                break
        r = f.readline()
        if all_seqs%1000==0:
            print "read %sM seqs" % (all_seqs/1e6)
    f.close()

    name_db = os.path.join(apa.path.polya_folder, "%s_pas" % poly_id)
    db.save(name_db+"_temp.bed.gz")

    f = gzip.open(name_db+"_temp.bed.gz")
    fout = gzip.open(name_db+"_temp.tab.gz", "wb")
    r = f.readline() # bed header
    r = f.readline()
    previous = None
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        chr, pos, cDNA = r[0], int(r[1]), int(r[-1])
        strand = "+" if cDNA>=0 else "-"
        pas = db_pas.get((chr, strand, pos), "")
        assert(len(pas)==1)
        distance = ""
        if previous!=None:
            if previous[0]==chr and previous[1]==strand:
                distance = pos - previous[-1]
        row = [r[0], r[1], r[3], list(pas)[0], distance]
        fout.write("\t".join(str(x) for x in row) + "\n")
        r = f.readline()
        previous = (chr, strand, pos)
    f.close()
    fout.close()

    print "filtering PAS database %s" % poly_id
    print "\tmin_distance = %snt" % min_distance
    print "\tPAS ranks = %s" % pas_signals

    def write_row(f, row):
        f.write("\t".join(str(el) for el in row) + "\n")

    f = gzip.open(name_db+"_temp.tab.gz")
    fout = gzip.open(name_db+"_temp_filtered.tab.gz", "wb")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        chr, pos, cDNA, pas, dist = r[0], int(r[1]), int(r[2]), r[3], r[4]
        dist = int(dist) if dist!="" else None
        if dist==None:
            pool = []
            write_row(fout, [chr, pos, cDNA, pas])
        elif dist>min_distance:
            # write pool?
            if len(pool)>1:
                pool = sorted(pool, key = lambda x : (x[0], -x[3])) # sort by pas rank (inc) + cDNA (dec)
                row = pool[0][1:-1]
                write_row(fout, row)
                pool = [(pas_signals.index(pas), chr, pos, cDNA, pas, dist)]
            elif len(pool)==0:
                pool = [(pas_signals.index(pas), chr, pos, cDNA, pas, dist)]
            elif len(pool)==1:
                row = pool[0][1:-1]
                write_row(fout, row)
                pool = [(pas_signals.index(pas), chr, pos, cDNA, pas, dist)]
        elif dist<min_distance:
            pool.append((pas_signals.index(pas), chr, pos, cDNA, pas, dist))
        r = f.readline()
    # write final record(s) from pool
    if len(pool)>0:
        pool.sort()
        write_row(fout, pool[0][1:-1])
    f.close()
    fout.close()

    b = pybio.data.Bedgraph()
    # add the distances to the new database
    f = gzip.open(name_db+"_temp_filtered.tab.gz")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        chr, pos, cDNA = r[0], int(r[1]), int(r[2])
        strand = "+" if cDNA>=0 else "-"
        b.set_value(chr, strand, pos, cDNA)
        r = f.readline()
    f.close()
    b.save(name_db+".bed.gz")

def read_polydb(poly_id):
    db = {}
    polyadb_tab = apa.path.polyadb_filename(poly_id, filetype="tab")
    if not os.path.exists(polyadb_tab):
        return db
    f = gzip.open(polyadb_tab)
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        data = dict(zip(header, r))
        chr = data["chr"]
        strand = data["strand"]
        pos = data["pos"]
        db["%s_%s_%s" % (chr, strand, pos)] = data
        r = f.readline()
    f.close()
    return db
