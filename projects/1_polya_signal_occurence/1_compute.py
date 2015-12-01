import apa
import pybio
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import numpy as np
import os
import cPickle

up = -60
down = 0

poly_id = "20150203_ina"
poly_id = "hg19"
genome = apa.polya.get_species(poly_id)
bg = pybio.data.Bedgraph(apa.path.polyadb_filename(poly_id, filetype="bed"))

def search_m(poly_id, motif_list = None):
    f = open("temp/%s.fasta" % poly_id, "wt")
    sites = 0
    vector = None
    for chr, strand_data in bg.raw.items():
        for strand, pos_data in strand_data.items():
            for pos, cDNA in pos_data.items():
                sites += 1
                seq = pybio.genomes.seq(genome, chr, strand, pos, start=up, stop=down)
                f.write(">%s\n%s\n" % (sites, seq))
                _, motif_vector = pybio.sequence.search(seq, motif_list)
                if len(motif_vector)==(down-up+1):
                    if vector==None:
                        vector = motif_vector
                    else:
                        vector = [x+y for x,y in zip(vector, motif_vector)]
    f.close()
    os.system("~/software/weblogo/seqlogo -f %s -o %s -F PNG -k 1 -e -c -w 70 -h 5 -Y -S -a -n -s -50" % ("temp/%s.fasta" % poly_id, "temp/%s" % poly_id))
    #return [x/float(sites) for x in vector], sites
    return vector, sites

def save_figure(data):
    mpl.rcParams['axes.labelsize'] = 8
    mpl.rcParams['axes.titlesize'] = 8
    mpl.rcParams['xtick.labelsize'] = 8
    mpl.rcParams['ytick.labelsize'] = 8
    mpl.rcParams['legend.fontsize'] = 8
    mpl.rc('axes', edgecolor='gray')
    mpl.rcParams['axes.linewidth'] = 0.1
    mpl.rcParams['legend.frameon'] = 'False'

    fig = plt.figure(figsize=(10, 10))
    axis = [plt.subplot(5,4,1)]
    for i in range(2, len(data)+1):
        axis.append(plt.subplot(5,4,i, sharex=axis[0]))

    #max_y = 0
    for aindex, (hexamer, vector) in enumerate(data):
        vector = pybio.utils.smooth(vector)
        axis[aindex].plot(range(0, len(vector)), vector, linewidth=1.3, alpha=1, antialiased=True)
        axis[aindex].xaxis.label.set_color("#a7a7a7")
        axis[aindex].yaxis.label.set_color("#a7a7a7")
        axis[aindex].get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        #max_y = max(max_y, max(vector))
        axis[aindex].legend([hexamer], loc='upper left')
    #plt.ylim(0, max_y)
    plt.xticks(range(0, len(vector), 10), range(up, down+1, 10))

    fig.tight_layout()
    fig.savefig("figures/polya_signal.png", dpi=300)
    fig.savefig("figures/polya_signal.pdf")

# main part

if not os.path.exists("temp"):
    os.makedirs("temp")
if not os.path.exists("figures"):
    os.makedirs("figures")

motif_list = ["AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATACA", "CATAAA", "AATATA", "GATAAA", "AATGAA", "AAGAAA", "ACTAAA", "AATAGA", "AATAAT", "AACAAA", "ATTACA", "ATTATA", "AACAAG", "AATAAG"]
#motif_list = ["AATAAA", "ATTAAA", "TATAAA", "AGTAAA"]

if not os.path.exists("data.pickle"):
    data = []
    legend = []
    for hexamer in motif_list:
        assert(len(hexamer)==6)
        vector, num_sites = search_m(poly_id, motif_list = [hexamer])
        print hexamer
        print vector
        print
        data.append((hexamer, vector))
        f = open("temp/%s.signal.tab" % hexamer, "wt")
        for freq in vector:
            f.write("%s\n" % freq)
        f.close()
    cPickle.dump(data, open("data.pickle", "wb"))
else:
    data = cPickle.load(open("data.pickle", "rb"))

save_figure(data)
