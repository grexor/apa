import pysam

bam_filename = "20120607_LUL24_e11_m1.bam"
bam_file = pysam.Samfile(bam_filename)
a_number = 0

ids = [67343995, 73450138, 117444518, 139702723, 141469021, 141961216]

for a in bam_file.fetch("2", 101622877, 101622886):
    
    # do not process spliced reads
    cigar = a.cigar
    cigar_types = [t for (t, v) in cigar]
    if 3 in cigar_types:
        continue
    
    read_id = int(a.qname)
    if read_id not in ids:
        continue
    
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
        
    print read_id, pos_end

