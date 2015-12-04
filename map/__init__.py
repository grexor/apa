import glob
import os
import apa
import pybio
import shutil
import commands

def map_experiment(lib_id, exp_id, map_id = 1, force=False):
    exp_id = int(exp_id)
    exp_data = apa.annotation.libs[lib_id].experiments[exp_id]
    map_folder = apa.path.map_folder(lib_id, exp_id, map_id = map_id)
    fastq_file = apa.path.map_fastq_file(lib_id, exp_id)
    #if exp_data["map_to"] not in ["hg19", "mm10", "dm5", "rn5"]:
    #    continue
    if os.path.exists(map_folder) and force==False:
        print "%s_e%s : MAP : skip (already mapped) or currently mapping" % (lib_id, exp_id)
        return
    if os.path.exists(map_folder):
        # remove folder
        if map_folder.find(lib_id)>0 and map_folder.find("m%s" % map_id)>0 and len(lib_id)>6: # security
            shutil.rmtree(map_folder)
    try:
        os.makedirs(map_folder)
    except:
        print "error creating mapping folder:", map_folder
        return False
    print "%s_e%s : MAP : %s" % (lib_id, exp_id, map_folder)
    pybio.map.star(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s" % (lib_id, exp_id, map_id))

def stats(lib_id):
    fname = os.path.join(apa.path.lib_folder(lib_id), "%s_stats.tab" % lib_id)
    print "writting to: %s" % fname
    f = open(fname, "wt")
    header = ["exp_id", "tissue", "condition", "replicate", "#reads [M]", "#mapped [M]" , "mapped [%]"]
    print "\t".join(header)
    f.write("\t".join(header) + "\n")
    for exp_id, exp_data in apa.annotation.libs[lib_id].experiments.items():
        fastq_file = apa.path.map_fastq_file(lib_id, exp_id)
        map_folder = apa.path.map_folder(lib_id, exp_id)
        bam_file = os.path.join(map_folder, "%s_e%s_m%s.bam" % (lib_id, exp_id, 1))
        # not fastq or bam perhaps? (bedgraph data)
        if not os.path.exists(fastq_file) or not os.path.exists(bam_file):
            continue
        num_reads = commands.getoutput("zcat %s | wc -l" % fastq_file).split("\n")[-1] # get last line of output
        num_reads = int(num_reads)/4
        map_reads = commands.getoutput("samtools view -c %s" % bam_file).split("\n")[-1] # get last line of output
        map_reads = int(map_reads)
        row = [exp_id, exp_data["tissue"], exp_data["condition"], exp_data["replicate"], "%.2f" % (num_reads/1e6), "%.2f" % (map_reads/1e6), "%.2f" % (map_reads*100.0/max(1, num_reads))]
        print "\t".join(str(x) for x in row)
        f.write("\t".join(str(x) for x in row) + "\n")
    f.close()
    return
