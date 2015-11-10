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

def search_m(poly_id, up=50, down=50, motif_list = None):
    genome = apa.polya.get_species(poly_id)
    f = open("temp/%s.fasta" % poly_id, "wt")
    print "processing: %s" % poly_id
    data = pybio.data.Bedgraph(apa.path.polyadb_filename(poly_id, filetype="bed"))
    sites = 0
    vector = [0]*(down+up+1)
    for chr, strand_data in data.raw.items():
        for strand, pos_data in strand_data.items():
            for pos, cDNA in pos_data.items():
                sites += 1
                seq = pybio.genomes.seq(genome, chr, strand, pos-up, pos+down)
                f.write(">%s\n%s\n" % (sites, seq))
                _, motif_vector = pybio.sequence.search(seq, motif_list)
                if len(motif_vector)==(down+up+1):
                    vector = [x+y for x,y in zip(vector, motif_vector)]
    f.close()
    os.system("~/software/weblogo/seqlogo -f %s -o %s -F PNG -k 1 -e -c -w 70 -h 5 -Y -S -a -n -s -50" % ("temp/%s.fasta" % poly_id, poly_id))
    print "%s, sites=%s" % (poly_id, sites)
    return [x/float(sites) for x in vector], sites

def save_figure(x, y, legend):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #ax.xaxis.label.set_color("#a7a7a7") # this changes axis ticks
    #ax.yaxis.label.set_color("#a7a7a7") # this changes axis ticks

    ax.tick_params(axis='x', colors='#a7a7a7') # this changes axis text
    ax.tick_params(axis='y', colors='#a7a7a7') # this changes axis text

    plt.xlim(x[0], x[-1])
    plt.ylim(0, 1)
    plt.xlabel("distance from poly-A site [nt]")
    plt.ylabel("signal ratio: present vs all")
    plt.title("Poly-A signal search")
    max_y = 0
    for vector in y:
        ax.plot(x, vector, linewidth=2, alpha=0.5, antialiased=True)
        max_y = max(max_y, max(vector))
    plt.ylim(0, max_y)
    plt.xticks([-50,-40,-30,-20,-10,0,10,20,30,40,50], [-50,-40,-30,-20,-10,0,10,20,30,40,50])

    ax.legend(legend, loc='upper right')
    fig.savefig("figures/polya_signal.png", dpi=300)
    fig.savefig("figures/polya_signal.pdf")

# main part

if not os.path.exists("temp"):
    os.makedirs("temp")
if not os.path.exists("figures"):
    os.makedirs("figures")

motif_list = ["AATAAA", "ATTAAA", "AGTAAA", "TATAAA"]
motif_list = ["AATAAA", "ATTAAA", "AGTAAA", "TATAAA"] # , "AAATAA", "ATAAAA", 'AAAATA', 'TTAAAA', 'TAAAAA', 'TAAAAT', 'ATAAAT', 'ATTAAA'] # 'AAAAAA'
motif_list = ["AATAAA", "ATAAAA", "AAATAA", "ATTAAA", "ATAAAT", "TAAAAA", "AAAATA", "TAATAA", "TAAAAT", "TTAAAA"] # NMF

x = range(-50, 51)
y = []
legend = []
for poly_id in ["20150311_miha", "hg19_tian", "hg19_derti", "20150203_ina", "hg19", "mm10", "rn5", "dm5"]:
    vector, num_sites = search_m(poly_id, motif_list = motif_list)
    y.append(vector)
    legend.append("%s (%s sites)" % (poly_id, '{0:,}'.format(num_sites)))
    f = open("temp/%s.signal.tab" % poly_id, "wt")
    for freq in vector:
        f.write("%s\n" % freq)
    f.close()

save_figure(x, y, legend)
