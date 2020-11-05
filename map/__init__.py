import glob
import os
import apa
import pybio
import sys
import json
import shutil
import regex
import gzip
from queue import *
from threading import *
from multiprocessing import Process # http://stackoverflow.com/questions/4496680/python-threads-all-executing-on-a-single-core

def map_experiment(lib_id, exp_id, map_id = 1, force=False, mapper="star", cpu=1, minlen=0.66, append=""):
    lib = apa.annotation.libs[lib_id]

    exp_id = int(exp_id)
    exp_data = apa.annotation.libs[lib_id].experiments[exp_id]
    map_folder = apa.path.map_folder(lib_id, exp_id, map_id = map_id, append=append)

    if os.path.exists(map_folder):
        if map_folder.find(lib_id)>0 and map_folder.find("m%s" % map_id)>0 and len(lib_id)>6: # security
            shutil.rmtree(map_folder)
    try:
        os.makedirs(map_folder)
    except:
        pass
    if os.path.exists(map_folder) and force==False:
        print("{lib_id}_e{exp_id} : MAP : skip (already mapped) or currently mapping".format(lib_id=lib_id, exp_id=exp_id))
        return

    if lib.seq_type=="single":
        fastq_file = apa.path.map_fastq_file(lib_id, exp_id, append=append)
        print("{lib_id}_e{exp_id} : MAP : {map_folder}".format(lib_id=lib_id, exp_id=exp_id, map_folder=map_folder))
        if mapper=="star":
            pybio.map.star(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s%s" % (lib_id, exp_id, map_id, append), cpu=cpu, minlen=minlen)
        if mapper=="sege":
            pybio.map.sege(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s%s" % (lib_id, exp_id, map_id, append), cpu=cpu)
        if mapper=="bowtie":
            pybio.map.bowtie(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s%s" % (lib_id, exp_id, map_id, append), cpu=cpu)
        if mapper=="bowtie2":
            pybio.map.bowtie2(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s%s" % (lib_id, exp_id, map_id, append), cpu=cpu)
        if mapper=="nano":
            pybio.map.nano(exp_data["map_to"], fastq_file, map_folder, "%s_e%s_m%s%s" % (lib_id, exp_id, map_id, append), cpu=cpu)

    if lib.seq_type=="paired":
        fastq_file1 = apa.path.map_fastq_file(lib_id, exp_id, append="_R1")
        fastq_file2 = apa.path.map_fastq_file(lib_id, exp_id, append="_R2")
        if mapper=="star":
            pybio.map.star_pair(exp_data["map_to"], fastq_file1, fastq_file2, map_folder, "%s_e%s_m%s%s" % (lib_id, exp_id, map_id, append), cpu=cpu, minlen=minlen)

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

def stats_to_tab(lib_id, map_id=1, append=""):
    fname_json = os.path.join(apa.path.lib_folder(lib_id), "%s_m%s%s.stats.json" % (lib_id, map_id, append))
    data = json.loads(open(fname_json).readline())
    fname = os.path.join(apa.path.lib_folder(lib_id), "%s_m%s%s.stats.tab" % (lib_id, map_id, append))
    print("writting statistics: {fname}".format(fname=fname))
    f = open(fname, "wt")
    header = [x for x in apa.annotation.libs[lib_id].experiments[list(apa.annotation.libs[lib_id].experiments.keys())[0]].keys() if x not in ["exp_id"]]
    header = ["exp_id"] + header + ["#reads [M]", "#mapped [M]" , "mapped [%]"]
    f.write("\t".join(header) + "\n")
    for exp_id, exp_data in apa.annotation.libs[lib_id].experiments.items():
        row = [exp_id]
        for x in header[1:-3]:
            row.append(exp_data[x])
        num_reads = data.get(str(exp_id), {}).get("num_reads", 0)
        map_reads = data.get(str(exp_id), {}).get("map_reads", 0)
        row = row + ["%.2f" % (num_reads/1e6), "%.2f" % (map_reads/1e6), "%.2f" % (map_reads*100.0/max(1, num_reads))]
        f.write("\t".join(str(x) for x in row) + "\n")
    f.close()

def stats(lib_id, map_id=1, append=""):
    for exp_id, exp_data in apa.annotation.libs[lib_id].experiments.items():
        stats_experiment(lib_id, exp_id, map_id=map_id)
    return

def stats_experiment(lib_id, exp_id, map_id=1, append=""):
    lib = apa.annotation.libs[lib_id]
    fname_json = os.path.join(apa.path.lib_folder(lib_id), "%s_m%s%s.stats.json" % (lib_id, map_id, append))
    map_folder = apa.path.map_folder(lib_id, exp_id, map_id=map_id, append=append)
    bam_file = os.path.join(map_folder, "%s_e%s_m%s%s.bam" % (lib_id, exp_id, map_id, append))
    if lib.seq_type=="single":
        fastq_files = [apa.path.map_fastq_file(lib_id, exp_id, append=append)]
    if lib.seq_type=="paired":
        fastq_files = [apa.path.map_fastq_file(lib_id, exp_id, append="_R1"), apa.path.map_fastq_file(lib_id, exp_id, append="_R2")]
    if os.path.exists(fname_json):
        data = json.loads(open(fname_json).readline())
    else:
        data = {}
    print("processing statistics: {lib_id}_e{exp_id}".format(lib_id=lib_id, exp_id=exp_id))
    if not os.path.exists(bam_file):
        return
    num_reads = 0
    for fastq_file in fastq_files:
        if not os.path.exists(fastq_file):
            continue
        output, error = pybio.utils.Cmd("bzcat {fastq_file} | wc -l".format(fastq_file=fastq_file)).run()
        temp_reads = output.split("\n")[0] # get last line of output
        temp_reads = int(temp_reads)/4
        num_reads += temp_reads
    output, error = pybio.utils.Cmd("samtools view -c {bam_file}".format(bam_file=bam_file)).run()
    map_reads = output.split("\n")[0] # get last line of output
    map_reads = int(map_reads)
    data[exp_id] = {"num_reads":num_reads, "map_reads":map_reads}
    f = open(fname_json, "wt")
    f.write(json.dumps(data))
    f.close()
    stats_to_tab(lib_id)
    return

# make transcript expression table (salmon)
def salmon(lib_id):
    library = apa.annotation.libs[lib_id]
    script_fname = os.path.join(apa.path.data_folder, lib_id, "salmon", "%s_salmon.sh" % lib_id)
    table_fname = os.path.join(apa.path.data_folder, lib_id, "%s_salmon.tab" % lib_id)
    salmon_folder = os.path.join(apa.path.data_folder, lib_id, "salmon")
    try:
        os.makedirs(salmon_folder)
    except:
        pass
    max_exp = len(library.experiments)
    f = open(script_fname, "wt")
    map_to = library.experiments[next(iter(library.experiments))]["map_to"]
    version = pybio.genomes.get_latest_version(map_to)
    genome_folder = os.path.join(pybio.path.genomes_folder, "%s.transcripts.%s.salmon" % (map_to, version))
    f.write("#!/bin/bash\n")
    for exp_id in library.experiments.keys():
        fastq_file = apa.path.map_fastq_file(lib_id, exp_id)
        output_folder = os.path.join(apa.path.data_folder, lib_id, "salmon", "e%s" % exp_id)
        f.write("salmon quant -i %s -l A -r <(bunzip2 -c %s) --validateMappings -o %s\n" % (genome_folder, fastq_file, output_folder))
    f.close()

    os.system("chmod +x %s" % script_fname)
    os.system(script_fname)

    # find name mapping between transcripts and genes
    transcript_gene = {}

    t_files = glob.glob(os.path.join(apa.path.pybio_folder, "genomes", "%s.transcripts.*/*.fa.gz" % library.genome))
    t_fname = t_files[0]
    f = pybio.data.Fasta(t_fname)

    print("reading transcript->gene data: {t_fname}".format(t_fname=t_fname))
    while f.read():
        id = f.id.split(" ")
        t_id = id[0]
        gene_id = id[3].split(":")[1]
        gene_name = ""
        if len(id)>6:
            gene_name = id[6].split(":")[1]
        transcript_gene[t_id] = [gene_id, gene_name]

    # combine results for individual experiments into the master quant table
    master_quant_fname = os.path.join(apa.path.data_folder, lib_id, "salmon", "%s_salmon.tab" % lib_id)
    fout = open(master_quant_fname, "wt")
    row = ["transcript_id", "gene_id", "gene_name"]
    for exp_id in library.experiments.keys():
        row.append("e%s_TPM" % exp_id)
    fout.write("\t".join(row) + "\n")
    data = {}
    for exp_id in library.experiments.keys():
        data[exp_id] = {}
        quant_fname = os.path.join(apa.path.data_folder, lib_id, "salmon", "e%s" % exp_id, "quant.sf")
        f = open(quant_fname)
        r = f.readline()
        r = f.readline()
        while r:
            r = r.replace("\n", "").replace("\r", "").split("\t")
            transcript_id = r[0]
            tpm = r[3]
            data[exp_id][transcript_id] = tpm
            r = f.readline()
        f.close()
    tids = set()
    for exp_id in library.experiments.keys():
        tids = tids.union(data[exp_id].keys())
    tids = list(tids)
    for tid in tids:
        gene_id, gene_name = transcript_gene[tid]
        row = [tid, gene_id, gene_name]
        for exp_id in library.experiments.keys():
            row.append(data[exp_id].get(tid, 0))
        fout.write("\t".join(row) + "\n")
    f.close()
