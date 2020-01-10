import glob
import pybio
import apa
import shutil
import os
import json
import datetime

files = glob.glob("%s/*/*.stats.tab" % apa.path.data_folder)
experiments = {}
reads = {}
methods = {}
genomes = set()

for fname in files:
    lib_id = fname.split("/")[-2]
    lib = apa.annotation.libs[lib_id]
    f = open(fname, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        data = dict(zip(header, r))
        nreads = float(data["#reads [M]"])
        exp_id = int(data["exp_id"])
        if "method" in data:
            method = data["method"]
        else:
            method = "not_speficied"
        try:
            genome = lib.experiments[exp_id]["map_to"]
            genomes.add(genome)
            experiments[genome] = experiments.get(genome, 0) + 1
            reads[genome] = reads.get(genome, 0) + nreads
            methods[method] = methods.get(method, 0) + 1
        except:
            pass
        r = f.readline()
    f.close()

today = datetime.date.today()
fname = os.path.join(apa.path.data_folder, "stats.json")
f = open(fname, "wt")
f.write(json.dumps({"experiments":experiments, "reads":reads, "methods":methods, "date":str(today)}))
f.close()

print("wrote statistics to file:", fname)
print(experiments)
print(reads)
print(methods)

# sort genomes by number of reads
data = []
for g in genomes:
    data.append((reads[g], g))
data.sort(reverse=True)

sum_experiments = 0
sum_reads = 0
print("Experiments\tGenome\tReads [M]")
for _, g in data:
    print("%s\t%s\t%s" % (experiments[g], g, reads[g]))
    sum_experiments+= experiments[g]
    sum_reads += reads[g]
print("All experiments = %s" % (sum_experiments))
print("All reads = %s M" % (sum_reads))
