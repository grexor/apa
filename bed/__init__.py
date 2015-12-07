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

def internal_priming(seq):
    if seq.count("A")<10:
        return False
    else:
        return True

def match_pas(seq):
    for hexamer in PAS_hexamers:
        if seq.find(hexamer)!=-1:
            return True
    return False

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

def bed_raw(lib_id, exp_id, map_id=1, force=False):
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
        apa.bed.bed_raw_paseq(lib_id, exp_id, map_id=1, force=force)
    if exp_data["method"]=="paseqx":
        apa.bed.bed_raw_paseqx(lib_id, exp_id, map_id=1, force=force)
    if exp_data["method"]=="lexfwd":
        apa.bed.bed_raw_lexfwd(lib_id, exp_id, map_id=1, force=force)

def bed_raw_paseq(lib_id, exp_id, map_id, force=False):
    assert(apa.annotation.libs[lib_id].experiments[exp_id]["method"] in ["pAseq", "paseq"])
    # http://www.cgat.org/~andreas/documentation/pysam/api.html
    # Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.

    r_filename = apa.path.r_filename(lib_id, exp_id)
    t_filename = apa.path.t_filename(lib_id, exp_id)

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
    for a in bam_file.fetch():
        a_number += 1

        if a_number%10000==0:
            sys.stdout.write("\r%s_e%s_m%s : %sK reads processed : %s (pas count = %s)" % (lib_id, exp_id, map_id, a_number/1000, bam_filename))
            sys.stdout.flush()

        # do not process spliced reads
        cigar = a.cigar
        cigar_types = [t for (t, v) in cigar]
        if 3 in cigar_types:
            continue

        read_id = int(a.qname)
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
        rnd_code = apa.annotation.rndcode(lib_id, read_id)

        aremoved = 0
        if a.is_reverse:
            last_cigar = a.cigar[0]
        else:
            last_cigar = a.cigar[-1]
        if last_cigar[0]==4:
            aremoved = last_cigar[1]

        key = "%s:%s" % (chr, strand)

        # update T file
        if aremoved>=6:
            true_site = True
            downstream_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=1, stop=15)
            upstream_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=-36, stop=-1)

            if downstream_seq.startswith("AAAA") or downstream_seq[:10].count("A")>=5 or upstream_seq.endswith("AAAA") \
                or upstream_seq[-10:].count("A")>=5:
                true_site = False

            if match_pas(upstream_seq):
                true_site = True

            if true_site:
                temp = dataT.get(key, {})
                temp2 = temp.get(pos_end, set())
                temp2.add(rnd_code)
                temp[pos_end] = temp2
                dataT[key] = temp

        # update R file
        temp = dataR.get(key, {})
        temp2 = temp.get(pos_end, set())
        temp2.add(rnd_code)
        temp[pos_end] = temp2
        dataR[key] = temp

    write_bed(dataR, r_filename)
    write_bed(dataT, t_filename)

def bed_raw_paseqx(lib_id, exp_id, map_id, force=False):
    assert(apa.annotation.libs[lib_id].experiments[exp_id]["method"]=="paseqx")

    r_filename = apa.path.r_filename(lib_id, exp_id)
    t_filename = apa.path.t_filename(lib_id, exp_id)

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
    for a in bam_file.fetch():
        a_number += 1

        if a_number%10000==0:
            print "%s_e%s_m%s : %sK reads processed : %s" % (lib_id, exp_id, map_id, a_number/1000, bam_filename)

        cigar = a.cigar
        cigar_types = [t for (t, v) in cigar]

        #if 3 in cigar_types: # skip spliced reads
        #    continue

        aremoved = 0
        if a.is_reverse:
            last_cigar = a.cigar[-1]
        else:
            last_cigar = a.cigar[0]
        if last_cigar[0]==4:
            aremoved = last_cigar[1]

        read_id = a.qname
        chr = bam_file.getrname(a.tid)
        strand = "+" if not a.is_reverse else "-"
        # we use the reference positions of the aligned read (aend, pos)
        # relative positions are stored in qend, qstart
        if strand=="+":
            pos_end = a.positions[0]
        else:
            pos_end = a.positions[-1]
        # for paseqx, we turn strand
        strand = {"+":"-", "-":"+"}[strand]

        key = "%s:%s" % (chr, strand)

        # update T file
        true_site = True
        downstream_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=1, stop=15)
        upstream_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=-36, stop=-1)

        if downstream_seq.startswith("AAAAA") or downstream_seq[:10].count("A")>=6:   #or upstream_seq.endswith("AAAA") or upstream_seq[-10:].count("A")>=5:
            true_site = False

        if match_pas(upstream_seq):
            true_site = True

        if true_site:
            temp = dataT.get(key, {})
            temp2 = temp.get(pos_end, set())
            temp2.add(read_id)
            temp[pos_end] = temp2
            dataT[key] = temp

        # update R file
        temp = dataR.get(key, {})
        temp2 = temp.get(pos_end, set())
        temp2.add(read_id)
        temp[pos_end] = temp2
        dataR[key] = temp

    # write R file
    write_bed(dataR, r_filename)
    # write T file
    write_bed(dataT, t_filename)

def bed_raw_lexfwd(lib_id, exp_id, map_id, force=False):
    assert(apa.annotation.libs[lib_id].experiments[exp_id]["method"]=="lexfwd")
    # http://www.cgat.org/~andreas/documentation/pysam/api.html
    # Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.

    r_filename = apa.path.r_filename(lib_id, exp_id)
    t_filename = apa.path.t_filename(lib_id, exp_id)

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
    pas_count = 0
    len_dist = {}
    ft = {}

    def aremoved_key(aremoved):
        if aremoved==0:
            return 0
        elif 1<=aremoved<=10:
            return 10
        elif 11<=aremoved<=20:
            return 20
        elif 21<=aremoved<=30:
            return 30
        elif 31<=aremoved<=40:
            return 40
        elif 41<=aremoved<=50:
            return 50
        elif 51<=aremoved<=60:
            return 60
        elif 61<=aremoved<=70:
            return 70
        elif 71<=aremoved<=80:
            return 80
        elif 81<=aremoved<=90:
            return 90
        elif 91<=aremoved<=100:
            return 100
        elif 101<=aremoved<=110:
            return 110
        elif 111<=aremoved<=120:
            return 120
        elif 121<=aremoved<=130:
            return 130
        elif 131<=aremoved<=140:
            return 140
        elif 141<=aremoved<=150:
            return 150

    for a in bam_file.fetch():
        a_number += 1

        #if a_number>10000:
        #    break

        if a_number%10000==0:
            print "%s_e%s_m%s : %sK reads processed : %s" % (lib_id, exp_id, map_id, a_number/1000, bam_filename)

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

        aremoved = 0
        if a.is_reverse:
            last_cigar = a.cigar[0]
        else:
            last_cigar = a.cigar[-1]
        if last_cigar[0]==4:
            aremoved = last_cigar[1]

        key = "%s:%s" % (chr, strand)

        t1 = ft.get(aremoved_key(aremoved), {})
        t2 = t1.get(key, {})
        t3 = t2.get(pos_end, set())
        t3.add(read_id)
        t2[pos_end] = t3
        t1[key] = t2
        ft[aremoved_key(aremoved)] = t1

        # search for 15A and update T files
        read_seq = a.seq # sequence of read
        sur_seq = pybio.genomes.seq(genome, chr, strand, pos_end, start=-18)

        #hit_read = regex.search(r'A(?:AAAAAAAAAAAAAAA){s<=4}', read_seq)
        internal = internal_priming(sur_seq)

        if not internal:
            temp = dataT.get(key, {})
            temp2 = temp.get(pos_end, set())
            temp2.add(read_id)
            temp[pos_end] = temp2
            dataT[key] = temp

        # update R file
        temp = dataR.get(key, {})
        temp2 = temp.get(pos_end, set())
        temp2.add(read_id)
        temp[pos_end] = temp2
        dataR[key] = temp
        len_dist[aremoved_key(aremoved)] = len_dist.get(aremoved_key(aremoved), 0) + 1

    write_bed(dataR, r_filename)
    write_bed(dataT, t_filename)

    f = open(os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_clipped_dist.tab" % lib_id), "wt")
    f.write("shows how many of the mapped reads were 3' clipped and by how much [nt]\n")
    f.write("\t".join(["clipped_3", "#reads", "percentage_of_all_mapped"])+"\n")
    keys = len_dist.keys()
    keys.sort()
    for k in keys:
        row = [k, len_dist[k], "%.2f" % (len_dist[k]/float(a_number))]
        f.write("\t".join([str(x) for x in row])+ "\n")
    f.close()

    for aremoved, data in ft.items():
        fname = os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.T%s.bed" % (lib_id, exp_id, map_id, aremoved))
        write_bed(data, fname)
    return

def bed_expression(lib_id, exp_id, map_id=1, force=False, polyid=None):
    """
    :param force: overwrite existing bedGraph files if True
    :param map_id: which mapping to take; default 1
    :param polyid: which poly-A atlas to use; if omitted the experiment species global atlas is used

    Generates expression bedGraph files for lib_id and exp_id.

    The bedGraph files are stored in:

    .. code-block:: bash

        ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/${lib_id}_e${exp_id}_m${map_id}.E.bed # E=expression

    """

    exp_id = int(exp_id)
    exp_data = apa.annotation.libs[lib_id].experiments[exp_id]
    map_to = exp_data["map_to"]
    if exp_data["method"] in ["pAseq", "paseq"]:
        apa.bed.bed_expression_paseq(lib_id, exp_id=exp_id, map_id=1, map_to=map_to, force=force)
    if exp_data["method"]=="paseqx":
        apa.bed.bed_expression_paseqx(lib_id, exp_id=exp_id, map_id=1, map_to=map_to, polyid=polyid, force=force)
        #apa.bed.bed_expression_lexpas(lib_id, exp_id=exp_id, map_id=1, map_to=map_to, polyid=polyid, force=force)
    if exp_data["method"]=="lexfwd":
        apa.bed.bed_expression_lexfwd(lib_id, exp_id=exp_id, map_id=1, map_to=map_to, polyid=polyid, force=force)
        #apa.bed.bed_expression_lexpas(lib_id, exp_id=exp_id, map_id=1, map_to=map_to, polyid=polyid, force=force)

def bed_expression_paseq(lib_id, exp_id, map_id, map_to, force=False):
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    r_filename = apa.path.r_filename(lib_id, exp_id)
    e_filename = apa.path.e_filename(lib_id, exp_id)
    e_filename_ucsc = apa.path.e_filename(lib_id, exp_id, filetype="ucsc")
    polyadb_filename = apa.path.polyadb_filename(genome)


    if os.path.exists(e_filename) and not force:
        print "%s_e%s_m%s : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s : E BED file : start" % (lib_id, exp_id, map_id)
        open(e_filename, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=-100, stop=25)
        e.save(e_filename)

    if os.path.exists(e_filename_ucsc) and not force:
        print "%s_e%s_m%s_ucsc : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s_ucsc : E BED file : start" % (lib_id, exp_id, map_id)
        open(e_filename_ucsc, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=-100, stop=25)
        e.save(e_filename_ucsc, genome=map_to, track_id="%s_e%s_m1" % (lib_id, exp_id))

def bed_expression_paseqx(lib_id, exp_id, map_id, map_to, polyid, force=False):
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    r_filename = apa.path.r_filename(lib_id, exp_id)

    if polyid==None:
        polyid = map_to
    polyadb_filename = apa.path.polyadb_filename(polyid)

    e_filename = apa.path.e_filename(lib_id, exp_id)
    if os.path.exists(e_filename) and not force:
        print "%s_e%s_m%s_ucsc : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s_ucsc : E BED file : start" % (lib_id, exp_id, map_id)
        open(e_filename, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=-100, stop=25)
        e.save(e_filename, track_id="%s_e%s_m1" % (lib_id, exp_id))

def bed_expression_lexfwd(lib_id, exp_id, map_id, map_to, polyid, force=False):
    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    r_filename = apa.path.r_filename(lib_id, exp_id)
    if polyid==None:
        polyid = map_to
    polyadb_filename = apa.path.polyadb_filename(polyid)

    e_filename = apa.path.e_filename(lib_id, exp_id)
    if os.path.exists(e_filename) and not force:
        print "%s_e%s_m%s_ucsc : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s_ucsc : E BED file : start" % (lib_id, exp_id, map_id)
        open(e_filename, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=-100, stop=25)
        e.save(e_filename, track_id="%s_e%s_m1" % (lib_id, exp_id))

def bed_expression_lexpas(lib_id, exp_id, map_id, map_to, polyid, force=False):

    region_start = 10
    region_stop = 60

    genome = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    r_filename = apa.path.r_filename(lib_id, exp_id)
    if polyid==None:
        polyid = map_to
    polyadb_filename = apa.path.polyadb_filename(polyid, filetype="pas")

    e_filename = apa.path.e_filename(lib_id, exp_id, filetype="pas")
    if os.path.exists(e_filename) and not force:
        print "%s_e%s_m%s_ucsc : E BED : already processed or currently processing" % (lib_id, exp_id, map_id)
    else:
        print "%s_e%s_m%s_ucsc : E BED file : start" % (lib_id, exp_id, map_id)
        open(e_filename, "wt").close() # touch E BED (processing)
        e = pybio.data.Bedgraph()
        e.overlay(polyadb_filename, r_filename, start=region_start, stop=region_stop)
        e.save(e_filename, track_id="%s_e%s_m1" % (lib_id, exp_id))
