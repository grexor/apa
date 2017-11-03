import pysam
import apa
import os
import pybio
import time
import glob
import sys
import regex

# http://www.cgat.org/~andreas/documentation/pysam/api.html
# Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.

PAS_hexamers = ['AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'CATAAA', 'GATAAA', 'AATATA', 'AATACA', 'AATAGA', 'ACTAAA', 'AAGAAA', 'AATGAA']

"""
25201104=D. Zheng and B. Tian, Systems Biology of RNA binding proteins. 2014
"This problem, commonly known as the "internal priming issue" can be partially
addressed by examining the genomic sequence surrounding the pA. For example,
our lab typically uses the -10 to +10 nt region, and considers a pA to be an
internal priming candidate if there are 6 continuous As or more than 7 As in
a 10 nt window within this region"
"""
def ip(s):
    if s.count("AAAAAA")>0:
        return True
    for i in range(0, len(s)):
        if s[i:i+10].count("A")>=7:
            return True
    return False

def match_pas(seq):
    for hexamer in PAS_hexamers:
        if seq.find(hexamer)!=-1:
            return True
    return False

def save(data, key, pos_end, read_id):
    level_1 = data.get(key, {})
    level_2 = level_1.get(pos_end, set())
    level_2.add(read_id)
    level_1[pos_end] = level_2
    data[key] = level_1

def write_bed(d, filename):
    """
    Save bedGraph file from dictionary d (chr->strand->position->value) to filename.
    """
    f = open(filename, "wt")
    for chr_strand, pos_data in d.items():
        chr, strand = chr_strand.split(":")
        positions = [(pos, len(rnd_set)) for pos, rnd_set in pos_data.items()]
        positions.sort()
        for pos, cDNA in positions:
            if strand=="+":
                f.write("%s\t%s\t%s\t%s\n" % (chr, pos, pos+1, cDNA))
            else:
                f.write("%s\t%s\t%s\t-%s\n" % (chr, pos, pos+1, cDNA))
    f.close()

def bed_raw(lib_id, exp_id, map_id=1, force=False, ip_filter=True):
    """
    :param force: overwrite existing bedGraph files if True
    :param map_id: which mapping to take; default 1
    Generates raw bedGraph files for lib_id and exp_id.

    The bedGraph files are stored in:

    .. code-block:: bash

        ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/${lib_id}_e${exp_id}_m${map_id}.R.bed # R=raw, unfiltered bedGraph files
        ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/${lib_id}_e${exp_id}_m${map_id}.T.bed # T=tail, filtered bedGraph files

    """
    lib = apa.annotation.libs[lib_id]
    exp_data = lib.experiments[exp_id]
    if exp_data["method"] in ["pAseq", "paseq"]:
        apa.bed.bed_raw_paseq(lib_id, exp_id, map_id=map_id, force=force, ip_filter=ip_filter)
    if exp_data["method"] in ["aseq"]:
        apa.bed.bed_raw_aseq(lib_id, exp_id, map_id=map_id, force=force, ip_filter=ip_filter)
    if exp_data["method"]=="lexrev":
        apa.bed.bed_raw_lexrev(lib_id, exp_id, map_id=map_id, force=force, ip_filter=ip_filter)
    if exp_data["method"]=="lexfwd":
        apa.bed.bed_raw_lexfwd(lib_id, exp_id, map_id=map_id, force=force, ip_filter=ip_filter)

def bed_raw_paseq(lib_id, exp_id, map_id, force=False, ip_filter=True):
    assert(apa.annotation.libs[lib_id].experiments[exp_id]["method"] in ["pAseq", "paseq"])
    # http://www.cgat.org/~andreas/documentation/pysam/api.html
    # Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.

    r_filename = apa.path.r_filename(lib_id, exp_id, map_id=map_id)
    t_filename = apa.path.t_filename(lib_id, exp_id, map_id=map_id)

    # don't redo analysis if files exists
    if (os.path.exists(r_filename) and not force) or (os.path.exists(t_filename) and not force):
        print "%s_e%s_m%s : R/T BED : already processed or currently processing" % (lib_id, exp_id, map_id)
        return

    lib = apa.annotation.libs[lib_id]
    exp_data = lib.experiments[exp_id]

    open(r_filename, "wt").close()
    open(t_filename, "wt").close()

    dataR = {}
    dataT = {}
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    bam_filename = os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.bam" % (lib_id, exp_id, map_id))
    bam_file = pysam.Samfile(bam_filename)
    a_number = 0
    ip_number = 0
    for a in bam_file.fetch():
        a_number += 1
        if a_number%10000==0:
            print "%s_e%s_m%s : %sM reads processed : %s" % (lib_id, exp_id, map_id, a_number/1e6, bam_filename)

        # do not process spliced reads
        #cigar = a.cigar
        #cigar_types = [t for (t, v) in cigar]
        #if 3 in cigar_types:
        #    continue

        read_id = int(a.qname)
        chr = bam_file.getrname(a.tid)
        strand = "+" if not a.is_reverse else "-"

        if strand=="+":
            pos_end = a.positions[-1] # aend points to one past the last aligned residue, also see a.positions
        else:
            pos_end = a.positions[0]
        rnd_code = apa.annotation.rndcode(lib_id, read_id)

        key = "%s:%s" % (chr, strand)

        save(dataR, key, pos_end, read_id)

        # internal priming
        check_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=-10, stop=10)
        if ip_filter:
            if ip(check_seq):
                ip_number += 1
                continue
        save(dataT, key, pos_end, read_id)

    write_bed(dataR, r_filename)
    write_bed(dataT, t_filename)

    f_ip = open(os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.ip_stats.txt" % (lib_id, exp_id, map_id)), "wt")
    f_ip.write("%s total processed alignments\n" % a_number)
    f_ip.write("%s (%.2f %%) alignments omitted due to internal priming\n" % (ip_number, ip_number/float(max(1, a_number))*100))
    f_ip.close()

def bed_raw_aseq(lib_id, exp_id, map_id, force=False, ip_filter=True):
    assert(apa.annotation.libs[lib_id].experiments[exp_id]["method"] in ["aseq"])
    # http://www.cgat.org/~andreas/documentation/pysam/api.html
    # Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.

    r_filename = apa.path.r_filename(lib_id, exp_id, map_id=map_id)
    t_filename = apa.path.t_filename(lib_id, exp_id, map_id=map_id)

    # don't redo analysis if files exists
    if (os.path.exists(r_filename) and not force) or (os.path.exists(t_filename) and not force):
        print "%s_e%s_m%s : R/T BED : already processed or currently processing" % (lib_id, exp_id, map_id)
        return

    lib = apa.annotation.libs[lib_id]
    exp_data = lib.experiments[exp_id]

    open(r_filename, "wt").close()
    open(t_filename, "wt").close()

    dataR = {}
    dataT = {}
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    bam_filename = os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.bam" % (lib_id, exp_id, map_id))
    bam_file = pysam.Samfile(bam_filename)
    a_number = 0
    ip_number = 0
    for a in bam_file.fetch():
        a_number += 1
        if a_number%10000==0:
            print "%s_e%s_m%s : %sM reads processed : %s" % (lib_id, exp_id, map_id, a_number/1e6, bam_filename)

        read_id = a.qname
        chr = bam_file.getrname(a.tid)
        strand = "+" if not a.is_reverse else "-"

        if strand=="+":
            pos_end = a.positions[-1] # aend points to one past the last aligned residue, also see a.positions
        else:
            pos_end = a.positions[0]

        key = "%s:%s" % (chr, strand)
        save(dataR, key, pos_end, read_id)
        check_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=-10, stop=10)
        if ip_filter:
            if ip(check_seq):
                ip_number += 1
                continue
        save(dataT, key, pos_end, read_id)

    write_bed(dataR, r_filename)
    write_bed(dataT, t_filename)

    f_ip = open(os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.ip_stats.txt" % (lib_id, exp_id, map_id)), "wt")
    f_ip.write("%s total processed alignments\n" % a_number)
    f_ip.write("%s (%.2f %%) alignments omitted due to internal priming\n" % (ip_number, ip_number/float(max(1, a_number))*100))
    f_ip.close()

def bed_raw_lexrev(lib_id, exp_id, map_id, force=False, ip_filter=True):
    assert(apa.annotation.libs[lib_id].experiments[exp_id]["method"]=="lexrev")
    read_len = apa.get_read_len(lib_id, exp_id)

    r_filename = apa.path.r_filename(lib_id, exp_id, map_id=map_id)
    t_filename = apa.path.t_filename(lib_id, exp_id, map_id=map_id)

    # don't redo analysis if files exists
    if (os.path.exists(t_filename) and not force):
        print "%s_e%s_m%s : tail BED : already processed or currently processing" % (lib_id, exp_id, map_id)
        return

    lib = apa.annotation.libs[lib_id]
    exp_data = lib.experiments[exp_id]

    open(t_filename, "wt").close()
    dataT = {}

    open(r_filename, "wt").close()
    dataR = {}

    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    bam_filename = os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.bam" % (lib_id, exp_id, map_id))
    bam_file = pysam.Samfile(bam_filename)
    a_number = 0
    ip_number = 0
    for a in bam_file.fetch():
        a_number += 1

        if a_number%10000==0:
            print "%s_e%s_m%s : %sM processed : %s" % (lib_id, exp_id, map_id, a_number/1e6, bam_filename)

        cigar = a.cigar
        cigar_types = [t for (t, v) in cigar]

        read_id = a.qname
        chr = bam_file.getrname(a.tid)
        strand = "+" if not a.is_reverse else "-"
        # we use the reference positions of the aligned read (aend, pos)
        # relative positions are stored in qend, qstart
        if strand=="+":
            pos_end = a.positions[0]
        else:
            pos_end = a.positions[-1]

        strand = {"+":"-", "-":"+"}[strand] # for lexrev, we turn strand
        key = "%s:%s" % (chr, strand)

        save(dataR, key, pos_end, read_id)

        # internal priming
        if ip_filter:
            check_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=-10, stop=10)
            if ip(check_seq):
                ip_number += 1
                continue

        save(dataT, key, pos_end, read_id)

    write_bed(dataR, r_filename)
    write_bed(dataT, t_filename)

    f_ip = open(os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.ip_stats.txt" % (lib_id, exp_id, map_id)), "wt")
    f_ip.write("%s total processed alignments\n" % a_number)
    f_ip.write("%s (%.2f %%) alignments omitted due to internal priming\n" % (ip_number, ip_number/float(max(1, a_number))*100))
    f_ip.close()

def bed_raw_lexfwd(lib_id, exp_id, map_id, force=False, ip_filter=True):
    assert(apa.annotation.libs[lib_id].experiments[exp_id]["method"]=="lexfwd")

    # http://www.cgat.org/~andreas/documentation/pysam/api.html
    # Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.

    t_filename = apa.path.t_filename(lib_id, exp_id, map_id=map_id)
    r_filename = apa.path.r_filename(lib_id, exp_id, map_id=map_id)

    # don't redo analysis if files exists
    if (os.path.exists(r_filename) and not force) or (os.path.exists(t_filename) and not force):
        print "%s_e%s_m%s : R/T BED : already processed" % (lib_id, exp_id, map_id)
        return

    lib = apa.annotation.libs[lib_id]
    exp_data = lib.experiments[exp_id]

    open(t_filename, "wt").close()
    open(r_filename, "wt").close()

    dataR = {}
    dataT = {}
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    bam_filename = os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.bam" % (lib_id, exp_id, map_id))
    bam_file = pysam.Samfile(bam_filename)

    a_number = 0
    ip_number = 0
    for a in bam_file.fetch():

        a_number += 1
        if a_number%100000==0:
            print "%s_e%s_m%s : %sM processed : %s, ip-filtering: %s" % (lib_id, exp_id, map_id, a_number/1e6, os.path.basename(bam_filename), ip_filter)

        read_id = a.qname
        chr = bam_file.getrname(a.tid)
        strand = "+" if not a.is_reverse else "-"
        # we use the reference positions of the aligned read (aend, pos)
        # relative positions are stored in qend, qstart
        if strand=="+":
            pos_end = a.aend - 1 # aend points to one past the last aligned residue, also see a.positions
            assert(pos_end==a.positions[-1])
        else:
            pos_end = a.pos
            assert(pos_end==a.positions[0])

        key = "%s:%s" % (chr, strand)

        save(dataR, key, pos_end, read_id)

        # internal priming
        if ip_filter:
            check_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=-10, stop=10)
            if ip(check_seq):
                ip_number += 1
                continue
        save(dataT, key, pos_end, read_id)

    write_bed(dataT, t_filename)
    write_bed(dataR, r_filename)

    f_ip = open(os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.ip_stats.txt" % (lib_id, exp_id, map_id)), "wt")
    f_ip.write("%s total processed alignments\n" % a_number)
    f_ip.write("%s (%.2f %%) alignments omitted due to internal priming\n" % (ip_number, ip_number/float(max(1, a_number))*100))
    f_ip.close()

    return

def bed_expression(lib_id, exp_id, map_id=1, force=False, poly_id=None, upstream=None, downstream=None):
    """
    :param force: overwrite existing bedGraph files if True
    :param map_id: which mapping to take; default 1
    :param poly_id: which poly-A atlas to use; if omitted the experiment species global atlas is used

    Generates expression bedGraph files for lib_id and exp_id.

    The bedGraph files are stored in:

    .. code-block:: bash

        ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/${lib_id}_e${exp_id}_m${map_id}.E.bed # E=expression

    """
    exp_id = int(exp_id)
    exp_data = apa.annotation.libs[lib_id].experiments[exp_id]
    map_to = exp_data["map_to"]

    upstream_defaults = {"pAseq":100, "paseq":100, "aseq":100, "lexrev":5, "lexfwd":100}
    downstream_defaults = {"pAseq":25, "paseq":25, "aseq":25, "lexrev":5, "lexfwd":25}
    if upstream==None:
        upstream = upstream_defaults[exp_data["method"]]
    if downstream==None:
        downstream = downstream_defaults[exp_data["method"]]

    if exp_data["method"] in ["pAseq", "paseq"]:
        apa.bed.bed_expression_paseq(lib_id, exp_id=exp_id, map_id=map_id, map_to=map_to, poly_id=poly_id, force=force, upstream=upstream, downstream=downstream)
    if exp_data["method"] in ["aseq"]:
        apa.bed.bed_expression_aseq(lib_id, exp_id=exp_id, map_id=map_id, map_to=map_to, poly_id=poly_id, force=force, upstream=upstream, downstream=downstream)
    if exp_data["method"]=="lexrev":
        apa.bed.bed_expression_lexrev(lib_id, exp_id=exp_id, map_id=map_id, map_to=map_to, poly_id=poly_id, force=force, upstream=upstream, downstream=downstream)
    if exp_data["method"]=="lexfwd":
        apa.bed.bed_expression_lexfwd(lib_id, exp_id=exp_id, map_id=map_id, map_to=map_to, poly_id=poly_id, force=force, upstream=upstream, downstream=downstream)

def bed_expression_paseq(lib_id, exp_id, map_id, map_to, poly_id, force=False, upstream=100, downstream=25):
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    r_filename = apa.path.r_filename(lib_id, exp_id, map_id=map_id)
    e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=poly_id)
    bam_filename = apa.path.bam_filename(lib_id, exp_id, map_id=map_id)
    if poly_id==None:
        poly_id = map_to
    polyadb_filename = apa.path.polyadb_filename(poly_id)

    e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=poly_id)
    if os.path.exists(e_filename) and not force:
        print "%s_e%s_m%s : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s : E BED, upstream=%s, downstream=%s" % (lib_id, exp_id, map_id, upstream, downstream)
        open(e_filename, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=-upstream, stop=downstream)
        #e.overlay2(polyadb_filename, bam_filename, start=-upstream, stop=downstream)
        e.save(e_filename, track_id="%s_e%s_m1" % (lib_id, exp_id))

def bed_expression_aseq(lib_id, exp_id, map_id, map_to, poly_id, force=False, upstream=100, downstream=25):
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    r_filename = apa.path.r_filename(lib_id, exp_id, map_id=map_id)
    e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=poly_id)
    bam_filename = apa.path.bam_filename(lib_id, exp_id, map_id=map_id)
    if poly_id==None:
        poly_id = map_to
    polyadb_filename = apa.path.polyadb_filename(poly_id)

    e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=poly_id)
    if os.path.exists(e_filename) and not force:
        print "%s_e%s_m%s : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s : E BED, upstream=%s, downstream=%s" % (lib_id, exp_id, map_id, upstream, downstream)
        open(e_filename, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=-upstream, stop=downstream)
        #e.overlay2(polyadb_filename, bam_filename, start=-upstream, stop=downstream)
        e.save(e_filename, track_id="%s_e%s_m1" % (lib_id, exp_id))

def bed_expression_lexrev(lib_id, exp_id, map_id, map_to, poly_id, force=False, upstream=5, downstream=5):
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    r_filename = apa.path.r_filename(lib_id, exp_id, map_id=map_id)
    if poly_id==None:
        poly_id = map_to
    polyadb_filename = apa.path.polyadb_filename(poly_id)

    bam_filename = apa.path.bam_filename(lib_id, exp_id, map_id=map_id)
    e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=poly_id)
    e_filename_norm = apa.path.e_filename_norm(lib_id, exp_id, map_id=map_id, poly_id=poly_id)
    if os.path.exists(e_filename) and not force:
        print "%s_e%s_m%s : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s : E BED, upstream=%s, downstream=%s" % (lib_id, exp_id, map_id, upstream, downstream)
        open(e_filename, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=-upstream, stop=downstream)
        #e.overlay2(polyadb_filename, bam_filename, start=-upstream, stop=downstream, reverse_strand=True)
        e.save(e_filename, track_id="%s_e%s_m1" % (lib_id, exp_id))
        e.norm()
        e.save(e_filename_norm, track_id="%s_e%s_m1" % (lib_id, exp_id), db_save="cpm")

def bed_expression_lexfwd(lib_id, exp_id, map_id, map_to, poly_id, force=False, upstream=100, downstream=25):
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    r_filename = apa.path.r_filename(lib_id, exp_id, map_id=map_id)
    if poly_id==None:
        poly_id = map_to
    polyadb_filename = apa.path.polyadb_filename(poly_id)

    bam_filename = apa.path.bam_filename(lib_id, exp_id, map_id=map_id)
    e_filename = apa.path.e_filename(lib_id, exp_id, map_id=map_id, poly_id=poly_id)
    if os.path.exists(e_filename) and not force:
        print "%s_e%s_m%s_ : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s : E BED, upstream=%s, downstream=%s" % (lib_id, exp_id, map_id, upstream, downstream)
        open(e_filename, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=-upstream, stop=downstream)
        #e.overlay2(polyadb_filename, bam_filename, start=-upstream, stop=downstream)
        e.save(e_filename, track_id="%s_e%s_m1" % (lib_id, exp_id))
