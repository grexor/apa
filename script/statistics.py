import glob
import pybio
import apa
import shutil
import os
import json

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

fname = os.path.join(apa.path.data_folder, "stats.json")
f = open(fname, "wt")
f.write(json.dumps({"experiments":experiments, "reads":reads, "methods":methods}))
f.close()

print "wrote statistics to file:", fname
print experiments
print reads
print methods

# sort genomes by number of reads
data = []
for g in genomes:
    data.append((reads[g], g))
data.sort(reverse=True)

sum_experiments = 0
sum_reads = 0
print "Experiments\tGenome\tReads [M]"
for _, g in data:
    print "%s\t%s\t%s" % (experiments[g], g, reads[g])
    sum_experiments+= experiments[g]
    sum_reads += reads[g]
print "All experiments = %s" % (sum_experiments)
print "All reads = %s M" % (sum_reads)

#import shutil
#fname_web = os.path.join(apa.config.expressrna_folder, "stats.json")
#shutil.copy(fname, fname_web)

"""
import matplotlib
matplotlib.use("Agg", warn=False)
import matplotlib.pyplot as plt
import math
import gzip
import copy
from matplotlib import cm as CM
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
c = mcolors.ColorConverter().to_rgb
matplotlib.rcParams['font.size'] = 17
matplotlib.rcParams['axes.labelsize'] = 17
matplotlib.rcParams['axes.titlesize'] = 17
matplotlib.rcParams['xtick.labelsize'] = 17
matplotlib.rcParams['ytick.labelsize'] = 17
matplotlib.rcParams['legend.fontsize'] = 17
matplotlib.rc('axes',edgecolor='gray')
matplotlib.rcParams['axes.linewidth'] = 0.1
matplotlib.rcParams['legend.frameon'] = 'False'

labels = ["%s (%s)" % (x,y) for x,y in zip(reads.keys(), experiments.values())]
sizes = experiments.values()
fig1, ax1 = plt.subplots()
cmap = plt.cm.Pastel2
colors = cmap(np.linspace(0., 1., len(reads)))
patches, texts, autotexts = ax1.pie(sizes, explode=[0]*len(reads), labels=labels, autopct='%1.1f%%', shadow=False, startangle=90, colors=colors)
for p in patches:
    p.set_edgecolor('white')
for autotext in texts:
    autotext.set_color('#555555')
for autotext in autotexts:
    autotext.set_color('#555555')
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
fig1.tight_layout()
fig1.savefig("1_statistics_experiments.png", bbox_inches='tight', dpi=300)

labels = ["%s (%s M)" % (x, int(y)) for x,y in zip(reads.keys(), reads.values())]
sizes = reads.values()
fig1, ax1 = plt.subplots()
cmap = plt.cm.Pastel2
colors = cmap(np.linspace(0., 1., len(reads)))
patches, texts, autotexts = ax1.pie(sizes, explode=[0]*len(reads), labels=labels, autopct='%1.1f%%', shadow=False, startangle=90, colors=colors)
for p in patches:
    p.set_edgecolor('white')
for autotext in texts:
    autotext.set_color('#555555')
for autotext in autotexts:
    autotext.set_color('#555555')
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
fig1.tight_layout()
fig1.savefig("1_statistics_reads.png", bbox_inches='tight', dpi=300)

shutil.copyfile("1_statistics_reads.png", "/home/gregor/expressrna/web/media/1_statistics_reads.png")
shutil.copyfile("1_statistics_experiments.png", "/home/gregor/expressrna/web/media/1_statistics_experiments.png")
"""
