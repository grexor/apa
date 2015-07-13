import pysam
import sys
import apa
import os
import pybio
import time
import glob
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['axes.titlesize'] = 10
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['legend.fontsize'] = 9

# all coordinates 0-based

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
            lib_id = "_".join(r.split("_")[:2])
            exp_id = int(r.split("_")[-1][1:])
            map_to = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
            species.add(map_to)
            r = f.readline()
        f.close()
        assert(len(species)==1)
        species = species.pop()
    return species

def process(poly_id):
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
            lib_id = "_".join(r.split("_")[:2])
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
        t_filename = apa.path.t_filename(lib_id, exp_id)
        if os.path.exists(t_filename):
            #bed.load(t_filename, meta=meta, min_cpm=min_cpm/float(10))
            bed.load(t_filename, meta=meta)
            print "%s: %s %s %s %s %.2fM" % (num_read, lib_id, exp_id, poly_id, os.path.exists(t_filename), bed.total_raw/1000000.0)

    # filter database
    bed.filter(min_distance=125)
    bed.save(apa.path.polyadb_filename(poly_id, filetype="temp"), db_save="raw")
    annotate(poly_id)
    #bed.save(apa.path.polyadb_filename(poly_id, filetype="bed"), db_save="raw", min_support=min_support)
    #bed.save(apa.path.polyadb_filename(poly_id, filetype="complete"), db_save="raw", min_support=min_support, filetype="complete")

def get_gene(species, gid):
    # only used by apa.annotate_position
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
    f = open(polyadb_temp, "rt")
    r = f.readline()
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        chr = r[0].replace("chr", "")
        pos = int(r[1])
        cDNA = float(r[-1])
        strand = "+" if cDNA>=0 else "-"
        gid_up, gid, gid_down, gid_interval = annotate_position(species, chr, strand, pos)
        key = "%s:%s:%s" % (chr, strand, pos)
        db[key] = (cDNA, gid, gid_interval)
        db_values.append(abs(cDNA))
        if gid!=None:
            gene_values.setdefault(gid, []).append(abs(cDNA))
        r = f.readline()
    f.close()

    cDNA_filter = {"hg19_derti":0, "hg19_tian":0}.get(poly_id, 10)
    accepted = {}
    for key, (cDNA, gid, gid_interval) in db.items():
        if abs(cDNA)>cDNA_filter:
            accepted[key] = 1

    ftab = open(polyadb_tab, "wt")
    ftab.write("\t".join(["chr", "strand", "pos", "gene_id", "gene_name", "interval", "cDNA", "seq_-125_0"]) + "\n")
    fbed = open(polyadb_bed, "wt")
    fbed.write("track type=bedGraph name=\"%s\" description=\"%s\" altColor=\"200,120,59\" color=\"120,101,172\" maxHeightPixels=\"100:50:0\" visibility=\"full\" priority=\"20\"" % (poly_id, poly_id))

    # store upstream sequences for polyar to classify weak/strong/other poly-A sites
    ffasta = open(polyadb_fasta, "wt")
    f = open(polyadb_temp, "rt")
    r = f.readline()
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        chr = r[0].replace("chr", "")
        pos = int(r[1])
        cDNA = float(r[-1])
        strand = "+" if cDNA>=0 else "-"
        # get upstream sequence
        if strand=="+":
            seq = pybio.genomes.seq(species, chr, strand, pos-125, pos)
        else:
            seq = pybio.genomes.seq(species, chr, strand, pos, pos+125)
        if seq==None:
            seq = "" # todo: check for which chrs there is no seq
        gid_up, gid, gid_down, gid_interval = annotate_position(species, chr, strand, pos)
        key = "%s:%s:%s" % (chr, strand, pos)
        if accepted.get(key, None)!=None:
            if gid==None:
                row = [chr, strand, pos, "", "", "", cDNA, seq]
            else:
                gene = get_gene(species, gid)
                interval = "%s:%s:%s" % (str(gid_interval[0]), str(gid_interval[1]), interval_types[gid_interval[2]])
                row = [chr, strand, pos, gid, gene["gene_name"], interval, cDNA, seq]
            ftab.write("\t".join(str(e) for e in row)+"\n")
            fbed.write("\t".join(r)+"\n")
            ffasta.write(">%s\n%s\n" % ("%s%s:%s" % (strand, chr, pos), seq))
        r = f.readline()
    f.close()
    ftab.close()
    fbed.close()
    ffasta.close()
    os.remove(polyadb_temp)

def annotate_pair(species, chr, strand, pos1, pos2):
    # keep positions in positive orientation
    if pos1>pos2:
        pos1, pos2 = pos2, pos1
    _, gid1, _, interval1 = annotate_position(species, chr, strand, pos1)
    _, gid2, _, interval2 = annotate_position(species, chr, strand, pos2)

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
    if itypes.find("io")==-1:
        return "tandem"
    if itypes.count("io")!=-1:
        # check upstream interval of upstream site
        if upint==None:
            return "skipped"
        if "%s%s" % (upint, itypes[0])=="oi":
            return "composite"
        return "skipped"
    return "other"

def annotate_position(species, chr, strand, pos, extension=5000):
    strand_os = "-" if strand=="+" else "+" # opposite strand
    gid_up, gid, gid_down, gid_interval = pybio.genomes.annotate(species, chr, strand, pos)
    default = (gid_up, gid, gid_down, gid_interval)
    if gid==None: # try to find gene, max upstream 5KB or middle of next gene on either strand
        new_pos = pos # in case checks fail, new position is old position
        gid_up_os, gid_os, gid_down_os, gid_interval_os = pybio.genomes.annotate(species, chr, strand_os, pos) # os = opposite strand
        if strand=="+":
            if gid_up==None or gid_os!=None: # no upstream gene OR gene on opposite strand present
                return default
            pos_up = get_gene(species, gid_up).get("gene_stop", None)
            pos_up_os = get_gene(species, gid_down_os).get("gene_stop", 0)
            if pos_up_os>pos_up: # there is a gene on the opposite strand before the upstream gene on the same strand?
                return default
            if abs(pos-pos_up)>extension:
                return default
            pos_down = get_gene(species, gid_down).get("gene_start", ())
            pos_down_os = get_gene(species, gid_up_os).get("gene_start", ())
            pos_down_min = min(pos_down, pos_down_os)
            if pos_down_min==(): # no downstream gene
                new_pos = get_gene(species, gid_up)["gene_intervals"][-1][1] # end of last interval
            else: # downstream gene
                if abs(pos_up-pos) <= abs(pos_down_min-pos_up)/2: # allow up to middle of downstream gene on either strand
                    new_pos = get_gene(species, gid_up)["gene_intervals"][-1][1] # end of last interval
        if strand=="-":
            if gid_up==None or gid_os!=None: # no upstream gene OR gene on opposite strand present
                return default
            pos_up = get_gene(species, gid_up).get("gene_start", None)
            pos_up_os = get_gene(species, gid_down_os).get("gene_start", ())
            if pos_up_os<pos_up: # there is a gene on the opposite strand before the upstream gene on the same strand?
                return default
            if abs(pos-pos_up)>extension:
                return default
            pos_down = get_gene(species, gid_down).get("gene_stop", 0)
            pos_down_os = get_gene(species, gid_up_os).get("gene_stop", 0)
            pos_down_max = max(pos_down, pos_down_os)
            if pos_down_max==0: # no downstream gene
                new_pos = get_gene(species, gid_up)["gene_intervals"][0][0] # start of first interval
            else: # downstream gene
                if abs(pos_up-pos) <= abs(pos_down_max-pos_up)/2: # allow up to middle of downstream gene on either strand
                    new_pos = get_gene(species, gid_up)["gene_intervals"][0][0] # start of first interval
        gid_up, gid, gid_down, gid_interval = pybio.genomes.annotate(species, chr, strand, new_pos)
    return gid_up, gid, gid_down, gid_interval
