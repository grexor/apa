import glob
import os
import apa
import pybio
import sys
import shutil
import commands
import regex
import gzip
from Queue import *
from threading import *
from multiprocessing import Process # http://stackoverflow.com/questions/4496680/python-threads-all-executing-on-a-single-core

def map_experiment(lib_id, exp_id, map_id = 1, force=False, mapper="star", cpu=1, minlen=0.66):
    exp_id = int(exp_id)
    exp_data = apa.annotation.libs[lib_id].experiments[exp_id]
    map_folder = apa.path.map_folder(lib_id, exp_id, map_id = map_id)
    fastq_file = apa.path.map_fastq_file(lib_id, exp_id)
    if os.path.exists(map_folder) and force==False:
        print "%s_e%s : MAP : skip (already mapped) or currently mapping" % (lib_id, exp_id)
        return
    if os.path.exists(map_folder):
        if map_folder.find(lib_id)>0 and map_folder.find("m%s" % map_id)>0 and len(lib_id)>6: # security
            shutil.rmtree(map_folder)
    try:
        os.makedirs(map_folder)
    except:
        pass

    print "%s_e%s : MAP : %s" % (lib_id, exp_id, map_folder)
    if mapper=="star":
        pybio.map.star(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s" % (lib_id, exp_id, map_id), cpu=cpu, minlen=minlen)
    if mapper=="sege":
        pybio.map.sege(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s" % (lib_id, exp_id, map_id), cpu=cpu)
    if mapper=="bowtie":
        pybio.map.bowtie(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s" % (lib_id, exp_id, map_id), cpu=cpu)

def preprocess_lexfwd(lib_id):
    threads = []
    for exp_id, exp_data in apa.annotation.libs[lib_id].experiments.items():
        t = Process(target=preprocess_lexfwd_thread, args=[exp_id, lib_id])
        threads.append(t)

    while len(threads)>0:
        to_process = min(apa.config.cores, len(threads))
        for i in range(to_process):
            threads[i].start()
        for i in range(to_process):
            threads[i].join()
        threads[:to_process] = []

def preprocess_lexfwd_thread(exp_id, lib_id):
    fastq_file_raw = apa.path.map_fastq_file_raw(lib_id, exp_id)
    fastq_file = apa.path.map_fastq_file(lib_id, exp_id)
    if not os.path.exists(fastq_file_raw):
        return
    processed, written, cut = 0, 0, 0
    fout = gzip.open(fastq_file, "wb")
    f = pybio.data.Fastq(fastq_file_raw)
    prog = regex.compile(r'AAA(?:AAAAAAA){s<=2}') # first AAA, then 7A with possible 2 mismatches
    while f.read():
        processed += 1
        i = prog.search(f.sequence)
        if i!=None:
            seq = f.sequence[:i.start(0)]
            cut += 1
        else:
            seq = f.sequence
        qual = f.quality[:len(seq)]
        if len(seq)>12:
            seq = seq[12:]
            qual = qual[12:]
        if len(seq)>10:
            fout.write("%s\n%s\n+\n%s\n" % (f.id, seq, qual))
            written += 1
        if processed%100000==0:
            sys.stdout.write("e%s, processed=%.2fM, cut=%s, kept=%s\n" % (exp_id, processed/1000000.0, "%.2f%%" % (float(cut)/processed*100.0), "%.2f%%" % (float(written)/processed*100.0)))
    fout.close()
    fname = os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "%s_e%s_preprocess.tab" % (lib_id, exp_id))
    f = open(fname, "wt")
    f.write("e%s\nprocessed=%.2fM\ncut=%s\nkept=%s" % (exp_id, processed/1000000.0, "%.2f%%" % (float(cut)/processed*100.0), "%.2f%%" % (float(written)/processed*100.0)))
    f.close()
    return

def stats(lib_id, map_id=1):
    fname = os.path.join(apa.path.lib_folder(lib_id), "%s_m%s.stats.tab" % (lib_id, map_id))
    print "writting to: %s" % fname
    f = open(fname, "wt")
    header = ["exp_id", "tissue", "condition", "replicate", "#reads [M]", "#mapped [M]" , "mapped [%]"]
    print "\t".join(header)
    f.write("\t".join(header) + "\n")
    for exp_id, exp_data in apa.annotation.libs[lib_id].experiments.items():
        fastq_file = apa.path.map_fastq_file(lib_id, exp_id)
        map_folder = apa.path.map_folder(lib_id, exp_id, map_id=map_id)
        bam_file = os.path.join(map_folder, "%s_e%s_m%s.bam" % (lib_id, exp_id, map_id))
        print fastq_file, bam_file
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
