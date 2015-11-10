import apa
import pybio
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import numpy as np
import os

mpl.rcParams['axes.labelsize'] = 15
mpl.rcParams['axes.titlesize'] = 12
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 12
mpl.rc('axes', edgecolor='gray')
mpl.rcParams['axes.linewidth'] = 0.3
mpl.rcParams['legend.frameon'] = 'False'

def save_figure(x, y, legend):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.tick_params(axis='x', colors='#a7a7a7') # this changes axis text
    ax.tick_params(axis='y', colors='#a7a7a7') # this changes axis text

    plt.xlim(x[0], x[-1])
    plt.xlabel("number of poly-A sites")
    plt.ylabel("number of genes")
    plt.title("poly-A site number per gene distribution")
    max_y = 0
    for vector in y:
        ax.plot(x, vector, linewidth=2, alpha=0.6, antialiased=True)
        max_y = max(max_y, max(vector))
    plt.xticks([1,2,3,4,5,6,7,8,9,10], [1,2,3,4,5,6,7,8,9,">10"])

    ax.legend(legend, loc='upper right')
    fig.savefig("polya_gene_count.png", dpi=300)
    fig.savefig("polya_gene_count.pdf")

x = range(1, 11)
y = []
legend = []

for poly_id in ["20150311_miha", "hg19_tian", "hg19_derti", "20150203_ina", "hg19", "mm10", "rn5", "dm5"]:
    print "processing", poly_id
    polyadb_ann_filename = apa.path.polyadb_filename(poly_id, filetype="tab")
    vector = [0] * 10
    f = open(polyadb_ann_filename, "rt")
    header = f.readline().replace("\n", "").replace("\r", "").split("\t")
    r = f.readline()
    gc = {}
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        data = dict(zip(header, r))
        if data["gene_id"]=="": # don't count not annotated positions
            r = f.readline()
            continue
        gc[data["gene_id"]] = gc.get(data["gene_id"], 0) + 1
        r = f.readline()
    f.close()

    for gid, count in gc.items():
        index = min(count-1, 9)
        vector[index] = vector[index] + 1

    y.append(vector)
    legend.append("%s" % (poly_id))
    f = open(os.path.join(apa.path.polya_folder, "%s.gene_count.tab" % poly_id), "wt")
    for freq in vector:
        f.write("%s\n" % freq)
    f.close()

save_figure(x, y, legend)
