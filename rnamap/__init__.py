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
import shutil

save_pdf = True

def read_deepbind(fname):
    rall = []
    rscaled = []
    f = open(fname, "rt")
    header = f.readline()
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        id = r[0][1:].split(" ")[0]
        gene_id, gene_name = id.split(":")

        # raw vectors
        v = []
        for el in r[1:]:
            v.append(float(el))
        while len(v)<400:
            v.append(0)
        rall.append(v)

        # min normalized vectors
        v = []
        for el in r[1:]:
            v.append(float(el))
        min_v = abs(min(v))
        while len(v)<400:
            v.append(-min_v)
        v = [e + min_v for e in v]
        v = [(gene_id, gene_name, sum(v))] + v
        rscaled.append(v)

        r = f.readline()
    f.close()

    # normalize
    r = []
    for nt in range(0, len(rall[0])):
        s = 0
        for i in range(0, len(rall)):
            s+=rall[i][nt]
        s = s / float(len(rall))
        r.append(s)
    return r[:-1], rscaled

def freq(data, all_genes):
    col_sums = []
    for i in range(1, 401+1):
        v = np.array(data[i]) # get column
        s = sum(v) #*0.9 # what vaue is 90% of the sum?
        col_sums.append(s)
    min_sum = 0.01 * max(col_sums) # don't show contribution for columns that have < 0.01 sum of the max
    temp = []
    for i in range(1, 401+1):
        v = np.array(data[i]) # get column
        v.sort() # sort
        v[:] = v[::-1] # reverse
        v = [x for x in v if x>0] # only consider >0 values otherwise cumsum will put cumulative values to 0 elements also
        s = sum(v) # *0.9 # what vaue is 90% of the sum?
        cs = np.cumsum(v) # cumulative sum
        res = [x for x in cs if x<=s] # how many elements (high->low) do i need to sum to get to 90% of the overall sum?
        if s>=min_sum and s>0:
            #print v
            #print cs
            #print res
            #print s
            #print len(res)
            #print
            #temp.append(len(res)/float(all_genes)) # append to frequency results
            temp.append(len(res)) # append #genes that explain 90% of the data at this position
        else:
            temp.append(0)
    temp = [e/float(all_genes)*100 for e in temp]
    return temp

def rnamap_deepbind(vpos, vneg, filename, title="test", ymax=None, site="proximal"):
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
    a = plt.axes([0.07, 0.2, 0.9, 0.7])
    a.set_xlim(0, 400)
    plt.ylabel("DeepBind score")
    plt.xlabel("distance (nt)")

    vpos_draw = pybio.utils.smooth(vpos)
    vneg_draw = pybio.utils.smooth(vneg)
    diff_draw = [x-y for x,y in zip(vpos_draw, vneg_draw)]
    plt.fill_between(range(0, len(diff_draw)), 0, diff_draw, facecolor='gray', alpha=0.6, interpolate=True)

    if ymax!=None:
        a.set_ylim(-ymax, ymax)

    p = mpatches.Rectangle([200, -200], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
    p = mpatches.Rectangle([0, 0], 400, ymax*0.001, facecolor='none', edgecolor=(0.8, 0.8, 0.8), linestyle='dotted')
    plt.gca().add_patch(p)
    plt.xticks([0,50,100,150,200,250,300,350,400], [-200,-150,-100,-50,0,50,100,150,200])

    print "saving", filename
    plt.title(title)
    plt.savefig(filename+".png", dpi=100)
    if save_pdf:
        plt.savefig(filename+".pdf")
    plt.close()

def rnamap_deepbind_heat(vpos, vneg, filename, title="test", site="proximal", stats=None, pair_type="tandem", alpha=0.8):
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
        if save_pdf:
            plt.savefig(filename+"_%s.pdf" % reg_type)
    return

def rnamap_area(vpos, vneg, vcon_up, vcon_down, filename, title="test", site="proximal", pair_type="tandem", ymax=None, stats=None):
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
    matplotlib.rcParams['axes.labelsize'] = 15
    matplotlib.rcParams['axes.titlesize'] = 15
    matplotlib.rcParams['xtick.labelsize'] = 13
    matplotlib.rcParams['ytick.labelsize'] = 13
    matplotlib.rcParams['legend.fontsize'] = 13
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.1
    matplotlib.rcParams['legend.frameon'] = 'False'
    from matplotlib import ticker

    fig = plt.figure(figsize=(20, 4))
    a = plt.axes([0.1, 0.2, 0.85, 0.7])
    a.set_xlim(0, 400)
    if ymax!=None:
        a.set_ylim(-ymax, ymax)
    plt.ylabel("counts per million (cpm)")
    plt.xlabel("distance [nt]")

    vpos = pybio.utils.smooth(vpos)
    vneg = [-el for el in pybio.utils.smooth(vneg)]
    vcon_up = pybio.utils.smooth(vcon_up)
    vcon_down = [-el for el in pybio.utils.smooth(vcon_down)]

    plt.fill_between(range(0, len(vpos)), 0, vpos, facecolor='red', alpha=0.6, interpolate=True)
    plt.plot(range(0, len(vpos)), vpos, color='red', alpha=1)

    plt.fill_between(range(0, len(vneg)), 0, vneg, facecolor='blue', alpha=0.6, interpolate=True)
    plt.plot(range(0, len(vneg)), vneg, color='blue', alpha=1)

    plt.plot(range(0, len(vcon_up)), vcon_up, color='black', alpha=1, linewidth=1.5) #, linestyle="dashed")
    plt.plot(range(0, len(vcon_down)), vcon_down, color='black', alpha=1, linewidth=1.5) #, linestyle="dashed")

    p = mpatches.Rectangle([200, -ymax], 0.01, ymax*2, facecolor='none', edgecolor=(0.8, 0, 0))
    plt.gca().add_patch(p)
    plt.xticks([0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400], [-200,-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200])
    plt.grid(alpha=0.2)
    print "saving", filename
    if site=="proximal":
        title = "%s.%s, e=%s, r=%s, c=%s" % (pair_type, site, stats["e.%s" % pair_type], stats["r.%s" % pair_type], stats["c_up.%s" % pair_type]+stats["c_down.%s" % pair_type])
    if site=="distal":
        title = "%s.%s, e=%s, r=%s, c=%s" % (pair_type, site, stats["r.%s" % pair_type], stats["e.%s" % pair_type], stats["c_up.%s" % pair_type]+stats["c_down.%s" % pair_type])
    plt.title(title)
    plt.savefig(filename+".png", dpi=100)
    if save_pdf:
        plt.savefig(filename+".pdf")
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
    matplotlib.rcParams['axes.labelsize'] = 17
    matplotlib.rcParams['axes.titlesize'] = 17
    matplotlib.rcParams['xtick.labelsize'] = 15
    matplotlib.rcParams['ytick.labelsize'] = 15
    matplotlib.rcParams['legend.fontsize'] = 15
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

        # A: divide by row sum
        #dn = d.ix[indices, 1:] # get only subset of rows
        #dn_div = dn.sum(axis=1) # get row sums
        #dn_div = dn_div.replace(0, 1) # do not divide row elements with 0
        #dn = dn.div(dn_div, axis=0) # divide by row sums

        # B: divide by max value in whole table
        dn = d.ix[indices, 1:]
        dn_div = dn.max()
        dn_div = dn_div.replace(0, 1)
        dn = dn.div(dn_div)

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
        labels = ["%s:%.2f" % (x[1], x[2]) for x in d.ix[indices, 0]]
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
        #plt.subplots_adjust(left=0.05, right=0.97, top=0.90, bottom=0.05) # nicely aligned version
        plt.subplots_adjust(left=0.1, right=0.97, top=0.90, bottom=0.05)
        cbar = fig.colorbar(heatmap, fraction=0.01, pad=0.01)
        print "saving %s" % (filename+"_%s.png" % reg_type)
        plt.savefig(filename+"_%s.png" % reg_type, dpi=150)
        if save_pdf:
            plt.savefig(filename+"_%s.pdf" % reg_type)
    return

def rnamap_freq(vpos, vneg, vcon_up, vcon_down, filename, title="test", site="proximal", stats=None, pair_type="tandem"):
    vpos = DataFrame(vpos)
    vneg = DataFrame(vneg)
    vcon_up = DataFrame(vcon_up)
    vcon_down = DataFrame(vcon_down)
    freq_pos = freq(vpos, len(vpos))
    freq_neg = freq(vneg, len(vneg))
    freq_con_up = freq(vcon_up, len(vcon_up))
    freq_con_down = freq(vcon_down, len(vcon_down))

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

    fig = plt.figure(figsize=(20, 4))
    a = plt.axes([0.1, 0.2, 0.9, 0.7])
    a.set_xlim(0, 400)
    plt.ylabel("%genes")
    plt.xlabel("distance [nt]")

    for axis in [a.xaxis, a.yaxis]:
        axis.set_major_locator(ticker.MaxNLocator(integer=True))

    vpos_graph = pybio.utils.smooth(freq_pos)
    vneg_graph = [-el for el in pybio.utils.smooth(freq_neg)]
    vcon_up_graph = pybio.utils.smooth(freq_con_up)
    vcon_down_graph = [-el for el in pybio.utils.smooth(freq_con_down)]

    ymax = max(max(vpos_graph), abs(min(vneg_graph)), max(vcon_up_graph), abs(min(vcon_down_graph)))
    ymax = math.ceil(ymax)

    a.set_ylim(-ymax, ymax)

    plt.fill_between(range(0, len(vpos_graph)), 0, vpos_graph, facecolor='red', alpha=0.1, interpolate=True)
    plt.plot(range(0, len(vpos_graph)), vpos_graph, color='red', alpha=1)

    plt.fill_between(range(0, len(vneg_graph)), 0, vneg_graph, facecolor='blue', alpha=0.1, interpolate=True)
    plt.plot(range(0, len(vneg_graph)), vneg_graph, color='blue', alpha=1)

    plt.plot(range(0, len(vcon_up_graph)), vcon_up_graph, color='black', alpha=1, linewidth=1.5) #, linestyle="dashed")
    plt.plot(range(0, len(vcon_down_graph)), vcon_down_graph, color='black', alpha=1, linewidth=1.5) #, linestyle="dashed")

    plt.margins(0)

    p = mpatches.Rectangle([200, -100], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
    plt.gca().add_patch(p)
    plt.xticks([0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400], [-200,-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200])
    if site=="proximal":
        title = "%s.%s, e=%s, r=%s, c=%s" % (pair_type, site, stats["e.%s" % pair_type], stats["r.%s" % pair_type], stats["c_up.%s" % pair_type]+stats["c_down.%s" % pair_type])
    if site=="distal":
        title = "%s.%s, e=%s, r=%s, c=%s" % (pair_type, site, stats["r.%s" % pair_type], stats["e.%s" % pair_type], stats["c_up.%s" % pair_type]+stats["c_down.%s" % pair_type])
    plt.title(title)
    plt.savefig(filename+".png", dpi=100)
    if save_pdf:
        plt.savefig(filename+".pdf")
    plt.close()

def process(comps_id=None, tab_file=None, clip_file="", genome=None, map_type="original", map_subfolder="rnamap", map_main_folder=None):

    if comps_id!=None:
        comps = apa.comps.read_comps(comps_id)
        genome = comps.species
        if comps.iCLIP_filename!=None:
            clip_file = os.path.join(apa.path.iCLIP_folder, comps.iCLIP_filename)
        tab_file = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
        rnamap_dest = os.path.join(apa.path.comps_folder, comps_id, map_subfolder)
        if comps.polya_db!=None:
            polydb = apa.polya.read(comps.polya_db)
    else:
        if map_main_folder==None:
            print "specify map root folder"
            return
        rnamap_dest = os.path.join(map_main_folder, map_subfolder)
        comps = apa.comps.Comps()

    assert(len(rnamap_dest)>=10)
    if os.path.exists(rnamap_dest):
        shutil.rmtree(rnamap_dest)
    os.makedirs(rnamap_dest)

    #if comps.iCLIP_filename==None:
    #    return

    fasta_files = {}
    for site in ["proximal", "distal"]:
        for pair_type in ["tandem", "composite", "skipped"]:
            for reg in ["r", "e"]:
                k = "%s_%s_%s" % (site, pair_type, reg)
                fname = os.path.join(rnamap_dest, k+".fasta")
                fasta_files[k] = open(fname, "wt")

    pc_thr = 0.1
    fisher_thr = 0.1
    pair_dist = 450

    if clip_file!="":
        clip = pybio.data.Bedgraph(clip_file)
    else:
        clip = None

    # r = repressed, e = enhanced, c = control
    stats = Counter()
    cdata_vectors = {} # individual vectors for heatmap
    sdata_vectors = {} # individual vectors for heatmap
    pdata_vectors = {} # individual vectors for heatmap
    cdata = {} # these are sums of vectors
    sdata = {} # these are sums of vectors
    pdata = {} # these are sums of vectors
    for pair_type in ["tandem", "composite", "skipped"]:
        for site in ["siteup", "sitedown"]:
            for reg in ["r", "e", "c_up", "c_down"]:
                cdata["%s.%s.%s" % (reg, pair_type, site)] = [0] * 401
                sdata["%s.%s.%s" % (reg, pair_type, site)] = [0] * 401
                pdata["%s.%s.%s" % (reg, pair_type, site)] = [0] * 401
                cdata_vectors["%s.%s.%s" % (reg, pair_type, site)] = []
                sdata_vectors["%s.%s.%s" % (reg, pair_type, site)] = []
                pdata_vectors["%s.%s.%s" % (reg, pair_type, site)] = []

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

        if abs(siteup_pos-sitedown_pos)<pair_dist:
            r = f.readline()
            continue

        if pc>0 and abs(pc)>pc_thr and fisher<fisher_thr:
            reg_siteup = "e"
        if pc<0 and abs(pc)>pc_thr and fisher<fisher_thr:
            reg_siteup = "r"
        if abs(pc)<pc_thr and fisher>fisher_thr:
            reg_siteup = "c_up" if pc>0 else "c_down"

        # also set reg_sitedown accordingly to reg_siteup
        reg_sitedown = {"e":"r", "r":"e", "c_up":"c_down", "c_down":"c_up", None:None}[reg_siteup]
        stats["%s.%s" % (reg_siteup, pair_type)] += 1

        if reg_siteup==None:
            r = f.readline()
            continue

        seq_up = pybio.genomes.seq(genome, chr, strand, siteup_pos, start=-200, stop=200)
        seq_down = pybio.genomes.seq(genome, chr, strand, sitedown_pos, start=-200, stop=200)

        if comps.polya_db!=None:
            if map_type in ["pas", "cs"]:
                siteup_pos += {"-":-1, "+":1}[strand] * int(polydb[(chr, strand, siteup_pos)]["%s_loci" % map_type])
                sitedown_pos += {"-":-1, "+":1}[strand] * int(polydb[(chr, strand, sitedown_pos)]["%s_loci" % map_type])
            if map_type in ["pas_manual"]:
                trim_seq_up = seq_up[:200].rfind("AATAAA")
                trim_seq_down = seq_down[:200].rfind("AATAAA")
                if trim_seq_up!=-1:
                    siteup_pos += {"-":1, "+":-1}[strand] * (200-trim_seq_up)
                if trim_seq_down!=-1:
                    sitedown_pos += {"-":1, "+":-1}[strand] * (200-trim_seq_down)

        if reg_siteup in ["e", "r"]:
            fasta_files["proximal_%s_%s" % (pair_type, reg_siteup)].write(">%s:%s %s%s:%s\n%s\n" % (gene_id, gene_name, strand, chr, siteup_pos, seq_up))
            fasta_files["distal_%s_%s" % (pair_type, reg_sitedown)].write(">%s:%s %s%s:%s\n%s\n" % (gene_id, gene_name, strand, chr, sitedown_pos, seq_down))

        # stringent: hw=25, hwt=10
        # TGTG
        for (rtype, seq, site) in [(reg_siteup, seq_up, "siteup"), (reg_sitedown, seq_down, "sitedown")]:
            _, z = pybio.sequence.search(seq, "TGTG")
            z = pybio.sequence.filter(z, hw=20, hwt=6)
            sdata["%s.%s.%s" % (rtype, pair_type, site)] = [x+y for x,y in zip(sdata["%s.%s.%s" % (rtype, pair_type, site)],z)]
            # [(gene_id, gene_name), 0, 0, 1, 1 ....]
            z_vector = [(gene_id, gene_name, sum(z))] + z
            sdata_vectors["%s.%s.%s" % (rtype, pair_type, site)].append(z_vector)
            assert(len(sdata["%s.%s.%s" % (rtype, pair_type, site)])==401)

        # AATAAA
        for (rtype, seq, site) in [(reg_siteup, seq_up, "siteup"), (reg_sitedown, seq_down, "sitedown")]:
            _, z = pybio.sequence.search(seq, "AATAAA")
            #z = pybio.sequence.filter(z)
            pdata["%s.%s.%s" % (rtype, pair_type, site)] = [x+y for x,y in zip(pdata["%s.%s.%s" % (rtype, pair_type, site)],z)]
            z_vector = [(gene_id, gene_name, sum(z))] + z
            pdata_vectors["%s.%s.%s" % (rtype, pair_type, site)].append(z_vector)
            assert(len(pdata["%s.%s.%s" % (rtype, pair_type, site)])==401)

        # CLIP
        if clip!=None:
            for (rtype, pos, site) in [(reg_siteup, siteup_pos, "siteup"), (reg_sitedown, sitedown_pos, "sitedown")]:
                z = []
                for index, x in enumerate(range(pos-200, pos+201)):
                    # all CLIP data is in UCSC genomic format of chromosome names?
                    z.append(clip.get_value("chr"+chr, strand, x, db="raw"))
                if strand=="-":
                    z.reverse()
                cdata["%s.%s.%s" % (rtype, pair_type, site)] = [x+y for x,y in zip(cdata["%s.%s.%s" % (rtype, pair_type, site)],z)]
                # [(gene_id, gene_name), 0, 0, 1, 1 ....]
                z_vector = [(gene_id, gene_name, sum(z))] + z
                cdata_vectors["%s.%s.%s" % (rtype, pair_type, site)].append(z_vector)
                assert(len(cdata["%s.%s.%s" % (rtype, pair_type, site)])==401)

        r = f.readline()
    f.close() # end of reading gene data

    print stats.items()

    for f in fasta_files.values():
        f.close()

    for pair_type in ["tandem", "composite", "skipped"]:
        n = max(1, stats["e.%s" % pair_type])
        cdata["e.%s.siteup" % pair_type] = [e/float(n) for e in cdata["e.%s.siteup" % pair_type]]
        cdata["r.%s.sitedown" % pair_type] = [e/float(n) for e in cdata["r.%s.sitedown" % pair_type]]
        sdata["e.%s.siteup" % pair_type] = [e/float(n) for e in sdata["e.%s.siteup" % pair_type]]
        sdata["r.%s.sitedown" % pair_type] = [e/float(n) for e in sdata["r.%s.sitedown" % pair_type]]
        pdata["e.%s.siteup" % pair_type] = [e/float(n) for e in pdata["e.%s.siteup" % pair_type]]
        pdata["r.%s.sitedown" % pair_type] = [e/float(n) for e in pdata["r.%s.sitedown" % pair_type]]
        n = max(1, stats["r.%s" % pair_type])
        cdata["r.%s.siteup" % pair_type] = [e/float(n) for e in cdata["r.%s.siteup" % pair_type]]
        cdata["e.%s.sitedown" % pair_type] = [e/float(n) for e in cdata["e.%s.sitedown" % pair_type]]
        sdata["r.%s.siteup" % pair_type] = [e/float(n) for e in sdata["r.%s.siteup" % pair_type]]
        sdata["e.%s.sitedown" % pair_type] = [e/float(n) for e in sdata["e.%s.sitedown" % pair_type]]
        pdata["r.%s.siteup" % pair_type] = [e/float(n) for e in pdata["r.%s.siteup" % pair_type]]
        pdata["e.%s.sitedown" % pair_type] = [e/float(n) for e in pdata["e.%s.sitedown" % pair_type]]

        n = max(1, stats["c_up.%s" % pair_type])
        cdata["c_up.%s.siteup" % pair_type] = [e/float(n) for e in cdata["c_up.%s.siteup" % pair_type]]
        cdata["c_down.%s.sitedown" % pair_type] = [e/float(n) for e in cdata["c_down.%s.sitedown" % pair_type]]
        n = max(1, stats["c_down.%s" % pair_type])
        cdata["c_down.%s.siteup" % pair_type] = [e/float(n) for e in cdata["c_down.%s.siteup" % pair_type]]
        cdata["c_up.%s.sitedown" % pair_type] = [e/float(n) for e in cdata["c_up.%s.sitedown" % pair_type]]

    cmax = {"tandem":0, "composite":0, "skipped":0}
    smax = {"tandem":0, "composite":0, "skipped":0}
    pmax = {"tandem":0, "composite":0, "skipped":0}
    for pair_type in ["tandem", "composite", "skipped"]:
        for reg in ["e", "r"]:
            for site in ["siteup", "sitedown"]:
                cmax[pair_type] = max(cmax[pair_type], max(cdata["%s.%s.%s" % (reg, pair_type, site)]))
                smax[pair_type] = max(smax[pair_type], max(sdata["%s.%s.%s" % (reg, pair_type, site)]))
                pmax[pair_type] = max(pmax[pair_type], max(pdata["%s.%s.%s" % (reg, pair_type, site)]))

    for pair_type in ["tandem", "composite", "skipped"]:
        for (site, site2) in [("siteup", "proximal"), ("sitedown", "distal")]:

            # clip
            if clip!=None:
                rnamap_area(cdata["e.%s.%s" % (pair_type, site)], cdata["r.%s.%s" % (pair_type, site)], cdata["c_up.%s.%s" % (pair_type, site)], cdata["c_down.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "clip.%s.%s" % (pair_type, site)), title="%s.%s" % (pair_type, site2), ymax=cmax[pair_type], site=site2, pair_type=pair_type, stats=stats)
                rnamap_freq(cdata_vectors["e.%s.%s" % (pair_type, site)], cdata_vectors["r.%s.%s" % (pair_type, site)], cdata_vectors["c_up.%s.%s" % (pair_type, site)], cdata_vectors["c_down.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "clip_freq.%s.%s" % (pair_type, site)), pair_type=pair_type, stats=stats, site=site2)
                rnamap_heat(cdata_vectors["e.%s.%s" % (pair_type, site)], cdata_vectors["r.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "clip_heat.%s.%s" % (pair_type, site)), pair_type=pair_type, stats=stats, site=site2, title=comps_id)

            # ug
            if "ug" in comps.rnamaps:
                rnamap_area(sdata["e.%s.%s" % (pair_type, site)], sdata["r.%s.%s" % (pair_type, site)], sdata["c_up.%s.%s" % (pair_type, site)], sdata["c_down.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "seq.%s.%s" % (pair_type, site)), title="%s.%s" % (pair_type, site2), ymax=smax[pair_type], site=site2, pair_type=pair_type, stats=stats)
                rnamap_freq(sdata_vectors["e.%s.%s" % (pair_type, site)], sdata_vectors["r.%s.%s" % (pair_type, site)], sdata_vectors["c_up.%s.%s" % (pair_type, site)], sdata_vectors["c_down.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "seq_freq.%s.%s" % (pair_type, site)), pair_type=pair_type, stats=stats, site=site2)
                rnamap_heat(sdata_vectors["e.%s.%s" % (pair_type, site)], sdata_vectors["r.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "seq_heat.%s.%s" % (pair_type, site)), pair_type=pair_type, stats=stats, site=site2, title=comps_id, alpha=0.3)

            # pas
            if "pas" in comps.rnamaps:
                rnamap_area(pdata["e.%s.%s" % (pair_type, site)], pdata["r.%s.%s" % (pair_type, site)], pdata["c_up.%s.%s" % (pair_type, site)], pdata["c_down.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "pas.%s.%s" % (pair_type, site)), title="%s.%s" % (pair_type, site2), ymax=pmax[pair_type], site=site2, pair_type=pair_type, stats=stats)
                rnamap_freq(pdata_vectors["e.%s.%s" % (pair_type, site)], pdata_vectors["r.%s.%s" % (pair_type, site)], pdata_vectors["c_up.%s.%s" % (pair_type, site)], pdata_vectors["c_down.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "pas_freq.%s.%s" % (pair_type, site)), pair_type=pair_type, stats=stats, site=site2)
                rnamap_heat(pdata_vectors["e.%s.%s" % (pair_type, site)], pdata_vectors["r.%s.%s" % (pair_type, site)], os.path.join(rnamap_dest, "pas_heat.%s.%s" % (pair_type, site)), pair_type=pair_type, stats=stats, site=site2, title=comps_id, alpha=0.3)

    # deep bind
    if comps.deepbind!=None:
        os.chdir(os.path.join(os.getenv("HOME"), "software/deepbind/"))
        for pair_type in ["tandem", "composite", "skipped"]:
            for site in ["proximal", "distal"]:
                for reg in ["r", "e"]:
                    sname = "%s_%s_%s.fasta" % (site, pair_type, reg)
                    dname = "%s_%s_%s_deepbind.tab" % (site, pair_type, reg)
                    sname = os.path.join(apa.path.comps_folder, comps_id, map_subfolder, sname)
                    dname = os.path.join(apa.path.comps_folder, comps_id, map_subfolder, dname)
                    os.system("./deepbind %s < \"%s\" > \"%s\"" % (comps.deepbind, sname, dname))
                    print "./deepbind %s < \"%s\" > \"%s\"" % (comps.deepbind, sname, dname)

        # ymax db
        ymax_db = {"tandem":0, "composite":0, "skipped":0}
        for pair_type in ["tandem", "composite", "skipped"]:
            for site in ["proximal", "distal"]:
                neg_name = "%s_%s_r_deepbind.tab" % (site, pair_type)
                neg_name = os.path.join(apa.path.comps_folder, comps_id, map_subfolder, neg_name)
                pos_name = "%s_%s_e_deepbind.tab" % (site, pair_type)
                pos_name = os.path.join(apa.path.comps_folder, comps_id, map_subfolder, pos_name)
                vneg, _ = read_deepbind(neg_name)
                vpos, _ = read_deepbind(pos_name)
                vneg = pybio.utils.smooth(vneg)
                vpos = pybio.utils.smooth(vpos)
                diff = [abs(x-y) for x,y in zip(vpos, vneg)]
                ymax_db[pair_type] = max(ymax_db[pair_type], max(diff))

        for pair_type in ["tandem", "composite", "skipped"]:
            for site in ["proximal", "distal"]:
                neg_name = "%s_%s_r_deepbind.tab" % (site, pair_type)
                neg_name = os.path.join(apa.path.comps_folder, comps_id, map_subfolder, neg_name)
                pos_name = "%s_%s_e_deepbind.tab" % (site, pair_type)
                pos_name = os.path.join(apa.path.comps_folder, comps_id, map_subfolder, pos_name)
                vneg, vneg_vectors = read_deepbind(neg_name)
                vpos, vpos_vectors = read_deepbind(pos_name)
                rnamap_deepbind(vpos, vneg, os.path.join(rnamap_dest, "%s_%s_deepbind" % (pair_type, site)), ymax=ymax_db[pair_type], site=site, title="%s %s %s (%s)" % (comps_id, site, pair_type, comps.deepbind))
                rnamap_deepbind_heat(vpos_vectors, vneg_vectors, os.path.join(rnamap_dest, "%s_%s_hdeepbind" % (pair_type, site)), site=site, title="%s %s %s (%s)" % (comps_id, site, pair_type, comps.deepbind))

    f = open(os.path.join(rnamap_dest, "index.html"), "wt")
    f.write("<html>\n")

    head = """

<head>
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
</head>

<script type="text/javascript" src="https://apa-db.org/software/highslide/highslide/highslide.js"></script>
<link rel="stylesheet" type="text/css" href="https://apa-db.org/software/highslide/highslide/highslide.css" />

<script type="text/javascript">
    hs.graphicsDir = 'https://apa-db.org/software/highslide/highslide/graphics/';
    hs.showCredits = false;
</script>

<script>

$(document).ready(function () {
  $(".slidingDiv").hide();

  $('.show_hide').click(function (e) {
    $(".slidingDiv").slideToggle("fast");
    var val = $(this).text() == "- details" ? "+ details" : "- details";
    $(this).hide().text(val).fadeIn("fast");
    e.preventDefault();
  });
});

</script>

<style>

.highslide img {
   border: 0px;
   outline: none;
}

a {
    text-decoration: none;
}

.show_hide {
  font-size: 13px;
}

</style>
"""

    f.write(head+"\n")
    f.write("<body>\n")

    body = """
    <a href="#" class="show_hide">+ details</a>
    """

    if comps_id==None:
        comparison = tab_file
        species = genome
    else:
        comparison = comps_id
        species = comps.species

    body += """<div style="font-size: 12px; padding-left: 10px;" class="slidingDiv">
    """ + comparison + """ : """ + species + """<br><br>
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
    """

    if clip!=None:
        body += "iCLIP = %s<br>" % os.path.basename(clip_file)

    body += """
    </div>
    <br>
    Control
    <div style="font-size: 12px; padding-left: 10px;">
    """

    def exp_print_out(L):
        body = ""
        for (comp_id, experiments, comp_name) in L:
            body += "%s<br>" % comp_id
            for exp_string in experiments:
                lib_id = "_".join(exp_string.split("_")[:-1]) # exp_string = "libid_e1"
                exp_id = int(exp_string.split("_")[-1][1:])
                condition = apa.annotation.libs[lib_id].experiments[exp_id]["condition"]
                method = apa.annotation.libs[lib_id].experiments[exp_id]["method"]
                replicate = apa.annotation.libs[lib_id].experiments[exp_id]["replicate"]
                body += "&nbsp;&nbsp;%s, %s, %s, rep=%s<br>" % (exp_string, condition, method, replicate)
        return body

    body += exp_print_out(comps.control)

    body += """
    </div>
    <br>
    Test
    <div style="font-size: 12px; padding-left: 10px;">
    """

    body += exp_print_out(comps.test)

    body += """
    </div>
    </div>
    <br>
    """
    f.write(body+"\n")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 12px;'>")

    if clip!=None:
        for t in ["tandem", "composite", "skipped"]:
            f.write("<tr><td align=center></td><td align=center>%s: proximal</td><td align=center>%s: distal</td></tr>\n" % (t, t))
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>iCLIP<br>enh=%s<br>rep=%s<br>con=%s</td>" % (t, stats["e.%s" % t], stats["r.%s" % t], stats["c_up.%s" % t]+stats["c_down.%s" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip.%s.siteup.png" % t, "clip.%s.siteup.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip.%s.sitedown.png" % t, "clip.%s.sitedown.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>iCLIP<br>targets</td>" % t)
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_freq.%s.siteup.png" % t, "clip_freq.%s.siteup.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_freq.%s.sitedown.png" % t, "clip_freq.%s.sitedown.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>iCLIP<br>top enh</td>" % t)
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_heat.%s.siteup_pos.png" % t, "clip_heat.%s.siteup_pos.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_heat.%s.sitedown_pos.png" % t, "clip_heat.%s.sitedown_pos.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>iCLIP<br>top rep</td>" % t)
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_heat.%s.siteup_neg.png" % t, "clip_heat.%s.siteup_neg.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("clip_heat.%s.sitedown_neg.png" % t, "clip_heat.%s.sitedown_neg.png" % t))
            f.write("</tr>")
            f.write("<tr><td><br><br></td><td><br><br></td><td><br><br></td></tr>")
            f.write("\n")

    if comps.deepbind!=None:
        for t in ["tandem", "composite", "skipped"]:
            f.write("<tr><td align=center></td><td align=center>%s: proximal</td><td align=center>%s: distal</td></tr>\n" % (t, t))
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>DeepBind<br><a href='http://tools.genes.toronto.edu/deepbind/%s/index.html' target=_db>%s</a><br>enh=%s<br>rep=%s</td>" % (t, comps.deepbind, comps.deepbind, stats["e.%s" % t], stats["r.%s" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("%s_proximal_deepbind.png" % t, "%s_proximal_deepbind.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("%s_distal_deepbind.png" % t, "%s_distal_deepbind.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>DeepBind<br><a href='http://tools.genes.toronto.edu/deepbind/%s/index.html' target=_db>%s</a><br>enh=%s<br>rep=%s</td>" % (t, comps.deepbind, comps.deepbind, stats["e.%s" % t], stats["r.%s" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("%s_proximal_hdeepbind_neg.png" % t, "%s_proximal_hdeepbind_neg.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("%s_distal_hdeepbind_neg.png" % t, "%s_distal_hdeepbind_neg.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>DeepBind<br><a href='http://tools.genes.toronto.edu/deepbind/%s/index.html' target=_db>%s</a><br>enh=%s<br>rep=%s</td>" % (t, comps.deepbind, comps.deepbind, stats["e.%s" % t], stats["r.%s" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("%s_proximal_hdeepbind_pos.png" % t, "%s_proximal_hdeepbind_pos.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("%s_distal_hdeepbind_pos.png" % t, "%s_distal_hdeepbind_pos.png" % t))
            f.write("</tr>")
            f.write("<tr><td><br><br></td><td><br><br></td><td><br><br></td></tr>")

    if "ug" in comps.rnamaps:
        for t in ["tandem", "composite", "skipped"]:
            f.write("<tr><td align=center></td><td align=center>%s: proximal</td><td align=center>%s: distal</td></tr>\n" % (t, t))
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>UG<br>enh=%s<br>rep=%s</td>" % (t, stats["e.%s" % t], stats["r.%s" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq.%s.siteup.png" % t, "seq.%s.siteup.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq.%s.sitedown.png" % t, "seq.%s.sitedown.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>UG<br>targets</td>" % t)
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_freq.%s.siteup.png" % t, "seq_freq.%s.siteup.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_freq.%s.sitedown.png" % t, "seq_freq.%s.sitedown.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>UG<br>top enh</td>" % t)
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_heat.%s.siteup_pos.png" % t, "seq_heat.%s.siteup_pos.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_heat.%s.sitedown_pos.png" % t, "seq_heat.%s.sitedown_pos.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>UG<br>top rep</td>" % t)
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_heat.%s.siteup_neg.png" % t, "seq_heat.%s.siteup_neg.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("seq_heat.%s.sitedown_neg.png" % t, "seq_heat.%s.sitedown_neg.png" % t))
            f.write("</tr>")
            f.write("<tr><td><br><br></td><td><br><br></td><td><br><br></td></tr>")
            f.write("\n")

    if "pas" in comps.rnamaps:
        for t in ["tandem", "composite", "skipped"]:
            f.write("<tr><td align=center></td><td align=center>%s: proximal</td><td align=center>%s: distal</td></tr>\n" % (t, t))
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>PAS<br>enh=%s<br>rep=%s</td>" % (t, stats["e.%s" % t], stats["r.%s" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("pas.%s.siteup.png" % t, "pas.%s.siteup.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("pas.%s.sitedown.png" % t, "pas.%s.sitedown.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>PAS<br>targets</td>" % t)
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("pas_freq.%s.siteup.png" % t, "pas_freq.%s.siteup.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("pas_freq.%s.sitedown.png" % t, "pas_freq.%s.sitedown.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>PAS<br>top enh</td>" % t)
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("pas_heat.%s.siteup_pos.png" % t, "pas_heat.%s.siteup_pos.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("pas_heat.%s.sitedown_pos.png" % t, "pas_heat.%s.sitedown_pos.png" % t))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>UG<br>top rep</td>")
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("pas_heat.%s.siteup_neg.png" % t, "pas_heat.%s.siteup_neg.png" % t))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("pas_heat.%s.sitedown_neg.png" % t, "pas_heat.%s.sitedown_neg.png" % t))
            f.write("</tr>")
            f.write("<tr><td><br><br></td><td><br><br></td><td><br><br></td></tr>")
            f.write("\n")

    f.write("\n")
    f.write("</table>")

    f.write("</body>")
    f.write("</html>\n")
    f.close()
