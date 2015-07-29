import matplotlib
matplotlib.use("Agg", warn=False)
import matplotlib.pyplot as plt
import math
import gzip
from matplotlib import cm as CM
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
c = mcolors.ColorConverter().to_rgb
import random
import pybio
import apa
import os
from collections import Counter
import cPickle
from pandas import DataFrame
import math
import matplotlib.patches as mpatches
import pickle

def freq(data, all_genes):
    col_sums = []
    for i in range(1, 401+1):
        v = np.array(data[i]) # get column
        s = sum(v)*0.9 # what vaue is 90% of the sum?
        col_sums.append(s)
    min_sum = 0.01 * max(col_sums) # don't show contribution for columns that have < 0.01 sum of the max
    temp = []
    for i in range(1, 401+1):
        v = np.array(data[i]) # get column
        v.sort() # sort
        v[:] = v[::-1] # reverse
        s = sum(v)*0.9 # what vaue is 90% of the sum?
        cs = np.cumsum(v) # cumulative sum
        res = [x for x in cs if x<=s] # how many elements (high->low) do i need to sum to get to 90% of the overall sum?
        if s>=min_sum and s>0:
            #temp.append(len(res)/float(all_genes)) # append to frequency results
            temp.append(len(res)) # append #genes that explain 90% of the data at this position
        else:
            temp.append(0)
    return temp

def rnamap_area(vpos, vneg, filename, title="test", ymax=None):
    """
    Draw RNA maps
    """
    import matplotlib
    matplotlib.use("Agg", warn=False)
    import matplotlib.pyplot as plt
    import math
    import gzip
    from matplotlib import cm as CM
    import matplotlib.patches as mpatches
    import matplotlib.ticker as mticker
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.colors as mcolors
    c = mcolors.ColorConverter().to_rgb

    # styling
    matplotlib.rcParams['axes.labelsize'] = 17
    matplotlib.rcParams['axes.titlesize'] = 17
    matplotlib.rcParams['xtick.labelsize'] = 14
    matplotlib.rcParams['ytick.labelsize'] = 14
    matplotlib.rcParams['legend.fontsize'] = 14
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.3
    matplotlib.rcParams['legend.frameon'] = 'False'

    fig = plt.figure(figsize=(20, 4))
    a = plt.axes([0.05, 0.2, 0.9, 0.7])
    a.set_xlim(0, 400)
    if ymax!=None:
        a.set_ylim(-ymax, ymax)
    plt.ylabel("cpm / #genes")
    plt.xlabel("distance from CS (nt)")

    vpos = pybio.utils.smooth(vpos)
    vneg = [-el for el in pybio.utils.smooth(vneg)]

    plt.fill_between(range(0, len(vpos)), 0, vpos, facecolor='red', alpha=0.6, interpolate=True)
    plt.plot(range(0, len(vpos)), vpos, color='red', alpha=1)

    plt.fill_between(range(0, len(vneg)), 0, vneg, facecolor='blue', alpha=0.6, interpolate=True)
    plt.plot(range(0, len(vneg)), vneg, color='blue', alpha=1)

    p = mpatches.Rectangle([200, -100], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
    plt.gca().add_patch(p)
    plt.xticks([0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400], [-200,-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200])

    print "saving", filename
    plt.title(title)
    plt.savefig(filename+".png", dpi=100)
    plt.savefig(filename+".svg")
    plt.close()

def rnamap_heat(vpos, vneg, filename, title="test", site="proximal", stats=None, pair_type="tandem", alpha=0.8):
    """
    Draw RNA heatmaps
    """
    def log_transform(x):
        if type(x)==numpy.float64:
            return math.log(x+1)
        else:
            return x

    # styling
    matplotlib.rcParams['axes.labelsize'] = 13
    matplotlib.rcParams['axes.titlesize'] = 13
    matplotlib.rcParams['xtick.labelsize'] = 11
    matplotlib.rcParams['ytick.labelsize'] = 11
    matplotlib.rcParams['legend.fontsize'] = 11
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.2
    matplotlib.rcParams['legend.frameon'] = 'False'
    import matplotlib.colors as mcolors

    def make_colormap(seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)

    c = mcolors.ColorConverter().to_rgb
    cred = make_colormap([c('white'), c('#ff8383'), 0.5, c('#ff8383'), c('red')])
    cblue = make_colormap([c('white'), c('#838cff'), 0.5, c('#838cff'), c('blue')])

    for reg_type, d in [("pos", vpos), ("neg", vneg)]:
        d = DataFrame(d)

        pickle.dump(d, open("heatmap.pickle", "wb"))

        order = d.sum(axis=1).order(ascending=False).index
        vals = [x[2] for x in d.ix[order, 0]]
        s = sum(vals)
        c = 0
        indices = []
        for i in order:
            z = d.ix[i, 0]
            c += z[2]
            indices.append(i)

        indices = indices[:20]

        # Plot it out
        fig, ax = plt.subplots()

        # -min, divide by (max-min)
        #dn = d.ix[indices, 1:].sub(d.ix[indices, 1:].min(axis=1), axis=0)
        #dn_div = dn.max(axis=1)-dn.min(axis=1)
        #dn_div = dn_div.replace(0, 1) # do not divide row elements with 0
        #dn = dn.div(dn_div, axis=0)

        # divide by row sum
        dn = d.ix[indices, 1:] # get only subset of rows
        dn_div = dn.sum(axis=1) # get row sums
        dn_div = dn_div.replace(0, 1) # do not divide row elements with 0
        dn = dn.div(dn_div, axis=0) # divide by row sums
        if reg_type=="pos":
            #heatmap = ax.pcolor(d.ix[indices, 1:], cmap=plt.cm.Reds, alpha=alpha)
            heatmap = ax.pcolor(dn, cmap=cred, alpha=alpha)
        else:
            #heatmap = ax.pcolor(d.ix[indices, 1:], cmap=plt.cm.Blues, alpha=alpha)
            heatmap = ax.pcolor(dn, cmap=cblue, alpha=alpha)

        fig = plt.gcf()
        fig.set_size_inches(30, 5)

        ax.set_yticks(np.arange(d.ix[indices, 1:].shape[0]) + 0.5, minor=False)
        #ax.set_frame_on(False)
        labels = ["%s : %.2f" % (x[1], x[2]) for x in d.ix[indices, 0]]
        ax.set_yticklabels(labels, minor=False)
        ax.set_xticklabels([-200,-150,-100,-50,0,50,100,150,200], minor=False)

        #ax.grid(False)
        ax.set_xlim(0, 401)
        ax.set_ylim(0, len(indices))
        ax.invert_yaxis()
        ax = plt.gca()
        for t in ax.xaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        for t in ax.yaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        p = mpatches.Rectangle([200, -100], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
        plt.gca().add_patch(p)
        plt.title("%s (top 20)" % title)
        plt.subplots_adjust(left=0.05, right=0.97, top=0.90, bottom=0.05)
        cbar = fig.colorbar(heatmap, fraction=0.01, pad=0.01)
        print "saving %s" % (filename+"_%s.png" % reg_type)
        plt.savefig(filename+"_%s.png" % reg_type, dpi=150)
        plt.savefig(filename+"_%s.svg" % reg_type)
    return

def rnamap_freq(vpos, vneg, filename, title="test", site="proximal", stats=None, pair_type="tandem"):
    c = mcolors.ColorConverter().to_rgb

    vpos = DataFrame(vpos)
    vneg = DataFrame(vneg)
    freq_pos = freq(vpos, len(vpos))
    freq_neg = freq(vneg, len(vneg))

    # styling
    matplotlib.rcParams['axes.labelsize'] = 11
    matplotlib.rcParams['axes.titlesize'] = 11
    matplotlib.rcParams['xtick.labelsize'] = 9
    matplotlib.rcParams['ytick.labelsize'] = 9
    matplotlib.rcParams['legend.fontsize'] = 9
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.1
    matplotlib.rcParams['legend.frameon'] = 'False'
    from matplotlib import ticker

    fig = plt.figure(figsize=(20, 2))
    a = plt.axes([0.05, 0.25, 0.9, 0.6])
    a.set_xlim(0, 400)
    plt.ylabel("#genes (90% cov.)")
    plt.xlabel("distance from poly-A site (nt)")

    for axis in [a.xaxis, a.yaxis]:
        axis.set_major_locator(ticker.MaxNLocator(integer=True))

    vpos_graph = pybio.utils.smooth(freq_pos)
    vneg_graph = [-el for el in pybio.utils.smooth(freq_neg)]

    ymax = max(max(vpos_graph), abs(min(vneg_graph)))
    #ymax = math.ceil(ymax*100)/100.0 # round up at the 2nd decimal
    #ymax += 0.01
    ymax = math.ceil(ymax)

    a.set_ylim(-ymax, ymax)

    plt.fill_between(range(0, len(vpos_graph)), 0, vpos_graph, facecolor='red', alpha=0.1, interpolate=True)
    plt.plot(range(0, len(vpos_graph)), vpos_graph, color='red', alpha=1)

    plt.fill_between(range(0, len(vneg_graph)), 0, vneg_graph, facecolor='blue', alpha=0.1, interpolate=True)
    plt.plot(range(0, len(vneg_graph)), vneg_graph, color='blue', alpha=1)

    plt.margins(0)

    p = mpatches.Rectangle([200, -100], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
    plt.gca().add_patch(p)
    plt.xticks([0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400], [-200,-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200])

    if site=="proximal":
        title = "%s.%s, enh=%s genes, sil=%s genes " % (pair_type, site, stats["e.%s" % pair_type], stats["r.%s" % pair_type])
    if site=="distal":
        title = "%s.%s, enh=%s genes, sil=%s genes " % (pair_type, site, stats["r.%s" % pair_type], stats["e.%s" % pair_type])
    plt.title(title)
    #plt.tight_layout() # doesnt work for this type of graphs
    plt.savefig(filename+".png", dpi=100)
    plt.savefig(filename+".svg")
    plt.close()

def process(comps_id=None, tab_file=None, clip_file="", genome=None, rnamap_dest="."):

    if comps_id!=None:
        comps = apa.comps.read_comps(comps_id)
        genome = comps.species
        if comps.iCLIP_filename!=None:
            clip_file = os.path.join(apa.path.iCLIP_folder, comps.iCLIP_filename)
        tab_file = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
        rnamap_dest = os.path.join(apa.path.comps_folder, comps_id, "rnamap")
        if comps.polya_db!=None:
            polydb = apa.polya.read(comps.polya_db)

    #if comps.iCLIP_filename==None:
    #    return

    pc_thr = 0.1
    fisher_thr = 0.1
    pair_dist = 450

    if clip_file!="":
        clip = pybio.data.Bedgraph(clip_file)
    else:
        clip = pybio.data.Bedgraph()

    # r = repressed, e = enhanced, c = control
    stats = Counter()
    cdata_vectors = {} # individual vectors for heatmap
    sdata_vectors = {} # individual vectors for heatmap
    cdata = {} # these are sums of vectors
    sdata = {} # these are sums of vectors
    for pair_type in ["tandem", "composite", "skipped"]:
        for site in ["siteup", "sitedown"]:
            for reg in ["r", "e"]:
                cdata["%s.%s.%s" % (reg, pair_type, site)] = [0] * 401
                sdata["%s.%s.%s" % (reg, pair_type, site)] = [0] * 401
                cdata_vectors["%s.%s.%s" % (reg, pair_type, site)] = []
                sdata_vectors["%s.%s.%s" % (reg, pair_type, site)] = []

    f = open(tab_file, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        chr = data["chr"]
        strand = data["strand"]
        gene_id = data["gene_id"]
        gene_name = data["gene_name"]
        siteup_pos = int(data["siteup_pos"])
        sitedown_pos = int(data["sitedown_pos"])
        pc = float(data["pc"])
        fisher = float(data["fisher"])
        pair_type = data["pair_type"]

        reg_siteup = None
        if pc>0 and abs(pc)>pc_thr and fisher<fisher_thr:
            reg_siteup = "e"
        if pc<0 and abs(pc)>pc_thr and fisher<fisher_thr:
            reg_siteup = "r"
        if reg_siteup in [None, "c"]:
            r = f.readline()
            continue
        if abs(siteup_pos-sitedown_pos)<pair_dist:
            r = f.readline()
            continue

        reg_sitedown = {"e":"r", "r":"e"}[reg_siteup]
        stats["%s.%s" % (reg_siteup, pair_type)] += 1

        seq_up = pybio.genomes.seq(genome, chr, strand, siteup_pos-200, siteup_pos+200)
        seq_down = pybio.genomes.seq(genome, chr, strand, sitedown_pos-200, sitedown_pos+200)

        # stringent: hw=25, hwt=10
        # UG proximal
        _, z = pybio.sequence.search(seq_up, "TGTG")
        z = pybio.sequence.filter(z, hw=20, hwt=6)
        sdata["%s.%s.siteup" % (reg_siteup, pair_type)] = [x+y for x,y in zip(sdata["%s.%s.siteup" % (reg_siteup, pair_type)],z)]
        # [(gene_id, gene_name), 0, 0, 1, 1 ....]
        z_vector = [(gene_id, gene_name, sum(z))] + z
        sdata_vectors["%s.%s.siteup" % (reg_siteup, pair_type)].append(z_vector)
        assert(len(sdata["%s.%s.siteup" % (reg_siteup, pair_type)])==401)

        # UG distal
        _, z = pybio.sequence.search(seq_down, "TGTG")
        z = pybio.sequence.filter(z, hw=20, hwt=6)
        sdata["%s.%s.sitedown" % (reg_sitedown, pair_type)] = [x+y for x,y in zip(sdata["%s.%s.sitedown" % (reg_sitedown, pair_type)],z)]
        # [(gene_id, gene_name), 0, 0, 1, 1 ....]
        z_vector = [(gene_id, gene_name, sum(z))] + z
        sdata_vectors["%s.%s.sitedown" % (reg_sitedown, pair_type)].append(z_vector)
        assert(len(sdata["%s.%s.sitedown" % (reg_sitedown, pair_type)])==401)

        # CLIP proximal
        z = []
        for index, x in enumerate(range(siteup_pos-200, siteup_pos+201)):
            z.append(clip.get_value("chr"+chr, strand, x, db="cpm"))
        if strand=="-":
            z.reverse()
        cdata["%s.%s.siteup" % (reg_siteup, pair_type)] = [x+y for x,y in zip(cdata["%s.%s.siteup" % (reg_siteup, pair_type)],z)]
        # [(gene_id, gene_name), 0, 0, 1, 1 ....]
        z_vector = [(gene_id, gene_name, sum(z))] + z
        cdata_vectors["%s.%s.siteup" % (reg_siteup, pair_type)].append(z_vector)
        assert(len(cdata["%s.%s.siteup" % (reg_siteup, pair_type)])==401)

        # CLIP distal
        z = []
        for index, x in enumerate(range(sitedown_pos-200, sitedown_pos+201)):
            z.append(clip.get_value("chr" + chr, strand, x, db="cpm"))
        if strand=="-":
            z.reverse()
        cdata["%s.%s.sitedown" % (reg_sitedown, pair_type)] = [x+y for x,y in zip(cdata["%s.%s.sitedown" % (reg_sitedown, pair_type)],z)]
        # [(gene_id, gene_name), 0, 0, 1, 1 ....]
        z_vector = [(gene_id, gene_name, sum(z))] + z
        cdata_vectors["%s.%s.sitedown" % (reg_sitedown, pair_type)].append(z_vector)
        assert(len(cdata["%s.%s.sitedown" % (reg_sitedown, pair_type)])==401)

        r = f.readline()
    f.close() # end of reading gene data

    if not os.path.exists(rnamap_dest):
        os.makedirs(rnamap_dest)

    for pair_type in ["tandem", "composite", "skipped"]:
        n = max(1, stats["e.%s" % pair_type])
        cdata["e.%s.siteup" % pair_type] = [e/float(n) for e in cdata["e.%s.siteup" % pair_type]]
        cdata["r.%s.sitedown" % pair_type] = [e/float(n) for e in cdata["r.%s.sitedown" % pair_type]]
        sdata["e.%s.siteup" % pair_type] = [e/float(n) for e in sdata["e.%s.siteup" % pair_type]]
        sdata["r.%s.sitedown" % pair_type] = [e/float(n) for e in sdata["r.%s.sitedown" % pair_type]]
        n = max(1, stats["r.%s" % pair_type])
        cdata["r.%s.siteup" % pair_type] = [e/float(n) for e in cdata["r.%s.siteup" % pair_type]]
        cdata["e.%s.sitedown" % pair_type] = [e/float(n) for e in cdata["e.%s.sitedown" % pair_type]]
        sdata["r.%s.siteup" % pair_type] = [e/float(n) for e in sdata["r.%s.siteup" % pair_type]]
        sdata["e.%s.sitedown" % pair_type] = [e/float(n) for e in sdata["e.%s.sitedown" % pair_type]]

    cmax = {"tandem":0, "composite":0, "skipped":0}
    smax = {"tandem":0, "composite":0, "skipped":0}
    for pair_type in ["tandem", "composite", "skipped"]:
        for reg in ["e", "r"]:
            cmax[pair_type] = max(cmax[pair_type], max(cdata["%s.%s.sitedown" % (reg, pair_type)]))
            cmax[pair_type] = max(cmax[pair_type], max(cdata["%s.%s.siteup" % (reg, pair_type)]))
            smax[pair_type] = max(smax[pair_type], max(sdata["%s.%s.sitedown" % (reg, pair_type)]))
            smax[pair_type] = max(smax[pair_type], max(sdata["%s.%s.siteup" % (reg, pair_type)]))
    for pair_type in ["tandem", "composite", "skipped"]:
        rnamap_area(sdata["e.%s.siteup" % pair_type], sdata["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "seq.%s.siteup" % pair_type), title="%s.proximal" % pair_type, ymax=smax[pair_type])
        rnamap_area(sdata["e.%s.sitedown" % pair_type], sdata["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "seq.%s.sitedown" % pair_type), title="%s.distal" % pair_type, ymax=smax[pair_type])
        rnamap_area(cdata["e.%s.siteup" % pair_type], cdata["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "clip.%s.siteup" % pair_type), title="%s.proximal" % pair_type, ymax=cmax[pair_type])
        rnamap_area(cdata["e.%s.sitedown" % pair_type], cdata["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "clip.%s.sitedown" % pair_type), title="%s.distal" % pair_type, ymax=cmax[pair_type])

        # freq
        rnamap_freq(cdata_vectors["e.%s.siteup" % pair_type], cdata_vectors["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "clip_freq.%s.siteup" % pair_type), pair_type=pair_type, stats=stats, site="proximal")
        rnamap_freq(cdata_vectors["e.%s.sitedown" % pair_type], cdata_vectors["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "clip_freq.%s.sitedown" % pair_type), pair_type=pair_type, stats=stats, site="distal")
        rnamap_freq(sdata_vectors["e.%s.siteup" % pair_type], sdata_vectors["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "seq_freq.%s.siteup" % pair_type), pair_type=pair_type, stats=stats, site="proximal")
        rnamap_freq(sdata_vectors["e.%s.sitedown" % pair_type], sdata_vectors["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "seq_freq.%s.sitedown" % pair_type), pair_type=pair_type, stats=stats, site="distal")

        # heat
        rnamap_heat(cdata_vectors["e.%s.siteup" % pair_type], cdata_vectors["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "clip_heat.%s.siteup" % pair_type), pair_type=pair_type, stats=stats, site="proximal", title=comps_id)
        rnamap_heat(cdata_vectors["e.%s.sitedown" % pair_type], cdata_vectors["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "clip_heat.%s.sitedown" % pair_type), pair_type=pair_type, stats=stats, site="distal", title=comps_id)
        rnamap_heat(sdata_vectors["e.%s.siteup" % pair_type], sdata_vectors["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "seq_heat.%s.siteup" % pair_type), pair_type=pair_type, stats=stats, site="proximal", title=comps_id, alpha=0.3)
        rnamap_heat(sdata_vectors["e.%s.sitedown" % pair_type], sdata_vectors["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "seq_heat.%s.sitedown" % pair_type), pair_type=pair_type, stats=stats, site="distal", title=comps_id, alpha=0.3)

    f = open(os.path.join(rnamap_dest, "index.html"), "wt")
    f.write("<html>\n")

    head = """

<head>
</head>

<script type="text/javascript" src="https://apa-db.org/software/highslide/highslide/highslide.js"></script>
<link rel="stylesheet" type="text/css" href="https://apa-db.org/software/highslide/highslide/highslide.css" />

<script type="text/javascript">
    hs.graphicsDir = 'https://apa-db.org/rnamotifs2/highslide/highslide/graphics/';
    hs.showCredits = false;
</script>

<style>

.highslide img {
   border: 0px;
   outline: none;
}

a {
    text-decoration: none;
}

</style>
"""

    f.write(head+"\n")

    f.write("<body>\n")

    body = """<div style="font-size: 12px;">
    """ + comps_id + """ : """ + comps.species + """<br><br>
    Description
    <div style="font-size: 12px; padding-left: 10px;">
    <font color=red>red = pc>0</font>
    <br>
    <font color=blue>blue = pc<0</font>
    <br>
    pc = [ control(proximal) / sum(control) ] - [ test(proximal) / sum(test) ]
    </div>
    <br>
    Parameters
    <div style="font-size: 12px; padding-left: 10px;">
    pc threshold = """ + str(pc_thr) + """
    <br>
    fisher threshold = """ + str(fisher_thr) + """
    <br>
    pair distance at least = """ + str(pair_dist) + """
    <br>
    iCLIP = """ + os.path.basename(clip_file) + """
    <br>
    </div>
    </div>
    <br>
    """
    f.write(body+"\n")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 12px;'>")

    for t in ["tandem", "composite", "skipped"]:
        f.write("<tr><td align=center></td><td align=center>proximal (%s)</td><td align=center>distal (%s)</td></tr>\n" % (t, t))
        f.write("<tr>")
        f.write("<td align=right valign=center>iCLIP<br><font color=red>e=%s</font><br><font color=blue>r=%s</font></td>" % (stats["e.tandem"], stats["r.tandem"]))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip.%s.siteup.png" % t, "clip.%s.siteup.png" % t))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip.%s.sitedown.png" % t, "clip.%s.sitedown.png" % t))
        f.write("</tr>")
        f.write("<tr>")
        f.write("<td align=right valign=center>iCLIP<br>cont.</td>")
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_freq.%s.siteup.png" % t, "clip_freq.%s.siteup.png" % t))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_freq.%s.sitedown.png" % t, "clip_freq.%s.sitedown.png" % t))
        f.write("</tr>")
        f.write("<tr>")
        f.write("<td align=right valign=center>iCLIP<br><font color=red>top 20 enh</font></td>")
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_heat.%s.siteup_pos.png" % t, "clip_heat.%s.siteup_pos.png" % t))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_heat.%s.sitedown_pos.png" % t, "clip_heat.%s.sitedown_pos.png" % t))
        f.write("</tr>")
        f.write("<tr>")
        f.write("<td align=right valign=center>iCLIP<br><font color=blue>top 20 rep</font></td>")
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_heat.%s.siteup_neg.png" % t, "clip_heat.%s.siteup_neg.png" % t))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_heat.%s.sitedown_neg.png" % t, "clip_heat.%s.sitedown_neg.png" % t))
        f.write("</tr>")
        f.write("<tr><td><br><br></td><td><br><br></td><td><br><br></td></tr>")
        f.write("\n")

    for t in ["tandem", "composite", "skipped"]:
        f.write("<tr><td align=center></td><td align=center>proximal (%s)</td><td align=center>distal (%s)</td></tr>\n" % (t, t))
        f.write("<tr>")
        f.write("<td align=right valign=center>UG<br><font color=red>e=%s</font><br><font color=blue>r=%s</font></td>" % (stats["e.tandem"], stats["r.tandem"]))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq.%s.siteup.png" % t, "seq.%s.siteup.png" % t))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq.%s.sitedown.png" % t, "seq.%s.sitedown.png" % t))
        f.write("</tr>")
        f.write("<tr>")
        f.write("<td align=right valign=center>UG<br>cont.</td>")
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_freq.%s.siteup.png" % t, "seq_freq.%s.siteup.png" % t))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_freq.%s.sitedown.png" % t, "seq_freq.%s.sitedown.png" % t))
        f.write("</tr>")
        f.write("<tr>")
        f.write("<td align=right valign=center>UG<br><font color=red>top 20 enh</font></td>")
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_heat.%s.siteup_pos.png" % t, "seq_heat.%s.siteup_pos.png" % t))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_heat.%s.sitedown_pos.png" % t, "seq_heat.%s.sitedown_pos.png" % t))
        f.write("</tr>")
        f.write("<tr>")
        f.write("<td align=right valign=center>UG<br><font color=blue>top 20 rep</font></td>")
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_heat.%s.siteup_neg.png" % t, "seq_heat.%s.siteup_neg.png" % t))
        f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_heat.%s.sitedown_neg.png" % t, "seq_heat.%s.sitedown_neg.png" % t))
        f.write("</tr>")
        f.write("<tr><td><br><br></td><td><br><br></td><td><br><br></td></tr>")
        f.write("\n")

    f.write("\n")
    f.write("</table>")

    f.write("</body>")
    f.write("</html>\n")
    f.close()
