import os
import sys
import apa
import pybio
import shutil
from collections import Counter

def process(comps_id):

    name_folder = "motifs"
    fasta_folder = os.path.join(apa.path.comps_folder, comps_id, name_folder, "fasta")

    if comps_id!=None:
        comps = apa.comps.read_comps(comps_id)
        genome = comps.species
        if comps.iCLIP_filename!=None:
            clip_file = os.path.join(apa.path.iCLIP_folder, comps.iCLIP_filename)
        tab_file = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
        dest_folder = os.path.join(apa.path.comps_folder, comps_id, name_folder)
        if comps.polya_db!=None:
            polydb = apa.polya.read(comps.polya_db)

    assert(len(dest_folder)>5)
    if os.path.exists(dest_folder):
        shutil.rmtree(dest_folder)
    os.makedirs(dest_folder)
    os.makedirs(fasta_folder)

    #if comps.iCLIP_filename==None:
    #    return

    fasta_files = {}
    for site in ["proximal", "distal"]:
        for pair_type in ["tandem", "composite", "skipped"]:
            for reg in ["r", "e", "c"]:
                k = "%s_%s_%s" % (site, pair_type, reg)
                fname = os.path.join(fasta_folder, k+".fasta")
                fasta_files[k] = open(fname, "wt")

    pc_thr = 0.1
    fisher_thr = 0.1
    pair_dist = 450
    control_pc_thr = 0.1
    control_fisher_thr = 0.2

    # r = repressed, e = enhanced, c = control
    stats = Counter()

    f = open(tab_file, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        chr = data["chr"]
        strand = data["strand"]
        gene_id = data["gene_id"]
        gene_name = data["gene_name"]
        siteup_pos = int(data["siteup_pos"])
        sitedown_pos = int(data["sitedown_pos"])

        pc = float(data["pc"])
        fisher = float(data["fisher"])
        pair_type = data["pair_type"]

        reg_siteup = None
        if pc>0 and abs(pc)>pc_thr and fisher<fisher_thr:
            reg_siteup = "e"
        elif pc<0 and abs(pc)>pc_thr and fisher<fisher_thr:
            reg_siteup = "r"
        elif abs(pc)<control_pc_thr and fisher>control_fisher_thr:
            reg_siteup = "c"

        if reg_siteup in [None]:
            r = f.readline()
            continue

        if abs(siteup_pos-sitedown_pos)<pair_dist:
            r = f.readline()
            continue

        reg_sitedown = {"e":"r", "r":"e", "c":"c"}[reg_siteup]
        stats["%s.%s" % (reg_siteup, pair_type)] += 1

        seq_up = pybio.genomes.seq(genome, chr, strand, siteup_pos, start=-200, stop=200)
        seq_down = pybio.genomes.seq(genome, chr, strand, sitedown_pos, start=-200, stop=200)

        fasta_files["proximal_%s_%s" % (pair_type, reg_siteup)].write(">%s:%s %s%s:%s\n%s\n" % (gene_id, gene_name, strand, chr, siteup_pos, seq_up))
        fasta_files["distal_%s_%s" % (pair_type, reg_sitedown)].write(">%s:%s %s%s:%s\n%s\n" % (gene_id, gene_name, strand, chr, sitedown_pos, seq_down))

        r = f.readline()
    f.close() # end of reading gene data

    for f in fasta_files.values():
        f.close()

    f_stats = open(os.path.join(dest_folder, "site_stats.tab"), "wt")
    f_stats.write("# r = repressed, e = enhanced, c = control\n")
    f_stats.write("# relative to proximal site; if proximal site = r, distal = e, and vice-versa; if proximal = c, distal is also c\n")
    for pair_type in ["tandem", "composite", "skipped"]:
        for reg in ["r", "e", "c"]:
            f_stats.write("%s\t%s\t%s\n" % (pair_type, reg, stats["%s.%s" % (reg, pair_type)]))
    f_stats.close()

    dreme(comps_id)


def dreme(comps_id):
    name_folder = "motifs"
    fasta_folder = os.path.join(apa.path.comps_folder, comps_id, name_folder, "fasta")
    dreme_folder = os.path.join(apa.path.comps_folder, comps_id, name_folder, "dreme")

    if os.path.exists(dreme_folder):
        shutil.rmtree(dreme_folder)
    os.makedirs(dreme_folder)

    for site in ["proximal", "distal"]:
        for pair_type in ["tandem", "composite", "skipped"]:
            fasta_pos = os.path.join(fasta_folder, "%s_%s_%s.fasta" % (site, pair_type, "e"))
            fasta_neg = os.path.join(fasta_folder, "%s_%s_%s.fasta" % (site, pair_type, "r"))
            fasta_con = os.path.join(fasta_folder, "%s_%s_%s.fasta" % (site, pair_type, "c"))
            dest_pos = os.path.join(dreme_folder, "%s_%s_e" % (site, pair_type))
            dest_neg = os.path.join(dreme_folder, "%s_%s_r" % (site, pair_type))
            os.system("dreme -p %s -n %s -maxk 6 -o %s -g 1000 -norc -e 0.5" % (fasta_pos, fasta_con, dest_pos))
            os.system("dreme -p %s -n %s -maxk 6 -o %s -g 1000 -norc -e 0.5" % (fasta_neg, fasta_con, dest_neg))

    #dreme -p proximal_tandem_e.fasta -n proximal_tandem_c.fasta -k 4 -o abc
