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
import bisect
import json

save_pdf = True

def freq(data, all_genes):
    if len(data)==0:
        return [0]
    col_sums = []
    for i in range(1, 401+1):
        v = np.array(data[i]) # get column
        s = sum(v) * 1 # what vaue is x% (x = 100% or 90% etc.) of the sum?
        col_sums.append(s)
    min_sum = 0.01 * max(col_sums) # don't show contribution for columns that have < 0.01 sum of the max
    temp = []
    for i in range(1, 401+1):
        v = np.array(data[i]) # get column
        v.sort()
        v[:] = v[::-1] # reverse
        v = [x for x in v if x>0] # only consider >0 values otherwise cumsum will put cumulative values to 0 elements also
        s = sum(v) # *0.9 # what vaue is 90% of the sum?
        cs = np.cumsum(v)
        res = [x for x in cs if x<=s] # how many elements (high->low) do i need to sum to get to 90% of the overall sum?
        if s>=min_sum and s>0:
            temp.append(len(res)) # append number of genes that explain x% of the data at this position
        else:
            temp.append(0)
    temp = [e/float(all_genes)*100 for e in temp] # show percentage of genes (/all_genes, * 100)
    return temp

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
    plt.xlabel("distance from polyA site [nt]")

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
        title = "%s.%s, e=%s, r=%s, c=%s" % (pair_type, site, stats[("enhanced", pair_type)], stats[("repressed", pair_type)], stats[("control_up", pair_type)]+stats[("control_down", pair_type)])
    if site=="distal":
        title = "%s.%s, e=%s, r=%s, c=%s" % (pair_type, site, stats[("repressed", pair_type)], stats[("enhanced", pair_type)], stats[("control_up", pair_type)]+stats[("control_down", pair_type)])
    plt.title(title)
    plt.savefig(filename+".png", dpi=100, transparent=True)
    if save_pdf:
        plt.savefig(filename+".pdf")
    plt.close()

    tab_data = {}
    tab_data["x"] = range(-200, 200)
    tab_data["vpos"] = list(vpos)
    tab_data["vneg"] = list(vneg)
    tab_data["cup"] = list(vcon_up)
    tab_data["cdown"] = list(vcon_down)
    tab_data["ymax"] = ymax
    tab_data["num_r"] = stats[("repressed", pair_type)]
    tab_data["num_e"] = stats[("enhanced", pair_type)]
    tab_data["num_cup"] = stats[("control_up", pair_type)]
    tab_data["num_cdown"] = stats[("control_down", pair_type)]

    f = open(filename+".tab", "wt")
    f.write(json.dumps(tab_data))
    f.close()

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
        if len(d)==0:
            continue
        d = DataFrame(d)
        #order = d.sum(axis=1).order(ascending=False).index
        order = d.sum(axis=1).sort_values(ascending=False).index
        vals = [x[2] for x in d.ix[order, 0]]
        s = sum(vals)
        c = 0

        fig, ax = plt.subplots()
        dn = d.ix[order[:20], 1:] # divide by max value per line
        dn_div = dn.max()
        dn_div = dn_div.replace(0, 1)
        dn = dn.div(dn_div)

        if reg_type=="pos":
            heatmap = ax.pcolor(dn, cmap=cred, alpha=alpha, vmin=0, vmax=1)
        else:
            heatmap = ax.pcolor(dn, cmap=cblue, alpha=alpha, vmin=0, vmax=1)

        fig = plt.gcf()
        fig.set_size_inches(30, 5)

        ax.set_yticks(np.arange(d.ix[order[:20], 1:].shape[0]) + 0.5, minor=False)
        #ax.set_frame_on(False)

        gene_list = []
        for x in d.ix[order, 0]:
            gene_id, gene_name, clip_sum = x[0], x[1], x[2]
            gene_list.append((gene_id, gene_name, clip_sum))

        labels = ["%s:%.2f" % (x[1], x[2]) for x in d.ix[order[:20], 0]]
        ax.set_yticklabels(labels, minor=False)
        ax.set_xticklabels([-200,-150,-100,-50,0,50,100,150,200], minor=False)

        #ax.grid(False)
        ax.set_xlim(0, 401)
        ax.set_ylim(0, len(order[:20]))
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
        plt.subplots_adjust(left=0.1, right=0.97, top=0.90, bottom=0.05)
        cbar = fig.colorbar(heatmap, fraction=0.01, pad=0.01)
        print "saving %s" % (filename+"_%s.png" % reg_type)
        plt.savefig(filename+"_%s.png" % reg_type, dpi=100, transparent=True)
        if save_pdf:
            plt.savefig(filename+"_%s.pdf" % reg_type)

        # write gene data for this image in the tab file
        f = open(filename+"_%s.tab" % reg_type, "wt")
        header = ["gene_id", "gene_name", "clip"]
        f.write("\t".join(header) + "\n")
        for (gene_id, gene_name, clip_sum) in gene_list:
            f.write("\t".join(str(el) for el in [gene_id, gene_name, clip_sum]) + "\n")
        f.close()

        tab_data = {}
        tab_data["data"] = dn.to_json(orient='values') # DataFrame to json
        tab_data["x"] = range(-200, 200)
        tab_data["ylabels"] = labels
        f = open(filename+"_%s_json.tab" % reg_type, "wt")
        f.write(json.dumps(tab_data))
        f.close()

    return

def save_control_binding(cup, cdown, filename, site="proximal", pair_type="tandem"):
    for reg_type, d in [("cup", cup), ("cdown", cdown)]:
        if len(d)==0:
            continue
        d = DataFrame(d)
        #order = d.sum(axis=1).order(ascending=False).index
        order = d.sum(axis=1).sort_values(ascending=False).index
        gene_list = []
        for x in d.ix[order, 0]:
            gene_id, gene_name, clip_sum = x[0], x[1], x[2]
            gene_list.append((gene_id, gene_name, clip_sum))

        # write gene data for this image in the tab file
        f = open(filename+"_%s.tab" % reg_type, "wt")
        header = ["gene_id", "gene_name", "clip"]
        f.write("\t".join(header) + "\n")
        for (gene_id, gene_name, clip_sum) in gene_list:
            f.write("\t".join(str(el) for el in [gene_id, gene_name, clip_sum]) + "\n")
        f.close()
    return

def rnamap_freq(vpos, vneg, vcon_up, vcon_down, filename=None, return_ymax=False, title="test", site="proximal", stats=None, pair_type="tandem", ymax=None):
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
    a = plt.axes([0.1, 0.2, 0.85, 0.7])
    a.set_xlim(0, 400)
    plt.ylabel("%genes")
    plt.xlabel("distance from polyA site [nt]")

    for axis in [a.xaxis, a.yaxis]:
        axis.set_major_locator(ticker.MaxNLocator(integer=True))

    vpos_graph = pybio.utils.smooth(freq_pos)
    vneg_graph = [-el for el in pybio.utils.smooth(freq_neg)]
    vcon_up_graph = pybio.utils.smooth(freq_con_up)
    vcon_down_graph = [-el for el in pybio.utils.smooth(freq_con_down)]

    if ymax==None:
        ymax = max(max(vpos_graph), abs(min(vneg_graph)), max(vcon_up_graph), abs(min(vcon_down_graph)))
        ymax = math.ceil(ymax)

    if return_ymax:
        return ymax

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
        title = "%s.%s, e=%s, r=%s, c=%s" % (pair_type, site, stats[("enhanced", pair_type)], stats[("repressed", pair_type)], stats[("control_up", pair_type)]+stats[("control_down", pair_type)])
    if site=="distal":
        title = "%s.%s, e=%s, r=%s, c=%s" % (pair_type, site, stats[("repressed", pair_type)], stats[("enhanced", pair_type)], stats[("control_up", pair_type)]+stats[("control_down", pair_type)])
    plt.title(title)
    plt.savefig(filename+".png", dpi=100, transparent=True)
    if save_pdf:
        plt.savefig(filename+".pdf")
    plt.close()

def coords(strand, proximal, distal, s1, s2, surr=200):
    """
    todo: write docs
    print coords("-", 1, 100, None, 1000, surr=100) -> [(1, 49, 100), (100, 100, 49), (None, 0, 0), (1000, 100, 100)]
    print coords("+", 1, 100, None, 1000, surr=100) -> [(1, 100, 49), (100, 49, 100), (None, 0, 0), (1000, 100, 100)]
    """
    def dist(L):
        R = []
        for x,y in zip(L, L[1:]):
            #R.append(y-x-1)
            R.append(y-x) # coords("-", 213864426, 213814546, 213864427, 213864427)
        return R
    sites = []
    indices = {}
    for site_pos in [proximal, distal, s1, s2]:
        if site_pos!=None:
            bisect.insort_left(sites, site_pos) # keep list of sites sorted
    sites_num = len(sites)
    distances = dist(sites)
    regions = {}
    regions["0_up"] = surr
    regions["%s_down" % (sites_num-1)] = surr
    for index, d in enumerate(distances):
        regions["%s_down" % index] = min(d/2, surr)
        regions["%s_up" % (index+1)] = min(d/2, surr)
    result = []
    for site_pos in [proximal, distal, s1, s2]:
        if site_pos!=None:
            i = sites.index(site_pos)
            if strand=="+":
                r = (site_pos, regions["%s_up" % i], regions["%s_down" % i])
            else:
                r = (site_pos, regions["%s_down" % i], regions["%s_up" % i])
        else:
            r = (None, 0, 0)
        result.append(r)
    return result

def adjust_len(vector, len_up, len_down, surr):
    # if string, return string
    # if list, return list
    if type(vector)==str:
        result_vector = list(vector)
    else:
        result_vector = vector
    result_vector = [0] * (surr-len_up) + result_vector + [0] * (surr-len_down)
    if type(vector)==str:
        return "".join([str(x) for x in result_vector])
    else:
        return result_vector

def presence_vector(vector, len_up, len_down, surr):
    # P = present, A = absent
    result = [0] * (surr-len_up) + [1]*len(vector) + [0] * (surr-len_down)
    assert(len(result)==401)
    return result

def process(comps_id, surr=200):

    comps = apa.comps.Comps(comps_id)
    genome = comps.species
    tab_file = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
    rnamap_dest = os.path.join(apa.path.comps_folder, comps_id, "rnamap")

    assert(len(rnamap_dest)>=10)
    if os.path.exists(rnamap_dest):
        shutil.rmtree(rnamap_dest)
    os.makedirs(rnamap_dest)

    # read clip data
    clip = {}
    for clip_name in comps.CLIP:
        clip_file = os.path.join(apa.path.iCLIP_folder, clip_name)
        clip[clip_name] = pybio.data.Bedgraph2(clip_file, fixed_cDNA=1) # set all cDNA values at all present positions to 1

    fasta_files = {}
    for site in ["proximal", "distal"]:
        for pair_type in ["tandem", "composite", "skipped", "all"]:
            for reg in ["repressed", "enhanced"]:
                k = "%s_%s_%s" % (site, pair_type, reg)
                fname = os.path.join(rnamap_dest, k+".fasta")
                fasta_files[(site, pair_type, reg)] = open(fname, "wt")

    # r = repressed, e = enhanced, c = control
    stats = Counter()
    stats_bysite = Counter()
    gene_list = {}
    cdata_vectors = {} # individual vectors for heatmap
    sdata_vectors = {}
    pdata_vectors = {}
    cdata = {} # these are sums of vectors
    sdata = {}
    pdata = {}
    present = {}
    adata = {} # all relevant data to store in a json file
    cgenes = {} # current gene count
    for pair_type in ["tandem", "composite", "skipped", "all"]:
        cgenes[pair_type] = 0
        adata[pair_type] = {}
        for site in ["proximal", "distal", "s1", "s2"]:
            for reg in ["repressed", "enhanced", "control_up", "control_down"]:
                for clip_name in comps.CLIP:
                    cdata[(clip_name, site, reg, pair_type)] = [0] * 401
                    cdata_vectors[(clip_name, site, reg, pair_type)] = []
                sdata[(site, reg, pair_type)] = [0] * 401
                pdata[(site, reg, pair_type)] = [0] * 401
                present[(site, reg, pair_type)] = [0] * 401
                sdata_vectors[(site, reg, pair_type)] = []
                pdata_vectors[(site, reg, pair_type)] = []

    pair_dist_stat = {}

    f = open(tab_file, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()

    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        chr = data["chr"]
        strand = data["strand"]
        gene_id = data["gene_id"]

        if (comps.exclusive_genes!=[]) and (gene_id not in comps.exclusive_genes):
            r = f.readline()
            continue

        if gene_id in comps.ignore_genes:
            r = f.readline()
            continue

        gene_name = data["gene_name"]
        proximal_pos = int(data["proximal_pos"])
        distal_pos = int(data["distal_pos"])
        pair_dist = abs(proximal_pos-distal_pos)
        if pair_dist<comps.pair_dist:
            r = f.readline()
            continue
        pair_dist_stat[pair_dist] = pair_dist_stat.get(pair_dist, 0) + 1

        if data["s1"]=="None":
            s1_pos = None
        else:
            s1_pos = int(data["s1"])
        if data["s2"]=="None":
            s2_pos = None
        else:
            s2_pos = int(data["s2"])

        pair_type = data["pair_type"]

        # gene selection is done with DEXSeq
        if comps.site_selection in ["DEX", "DEX2", "DEX3", "DEX4"]:
            proximal_reg = data["gene_class"]

        # also set reg_distal accordingly to reg_proximal
        distal_reg = {"enhanced":"repressed", "repressed":"enhanced", "control_up":"control_down", "control_down":"control_up", None:None}[proximal_reg]
        if pair_type in ["tandem", "composite"]:
            s1_reg = s2_reg = distal_reg
        elif pair_type=="skipped":
            s1_reg = proximal_reg
            s2_reg = distal_reg

        if proximal_reg==None:
            r = f.readline()
            continue

        # store adata
        adata[pair_type][cgenes[pair_type]] = {}
        adata[pair_type][cgenes[pair_type]]["proximal_reg"] = proximal_reg; adata[pair_type][cgenes[pair_type]]["distal_reg"] = distal_reg; adata[pair_type][cgenes[pair_type]]["s1_reg"] = s1_reg; adata[pair_type][cgenes[pair_type]]["s2_reg"] = s2_reg
        for site_type in ["proximal", "distal", "s1", "s2"]:
            adata[pair_type][cgenes[pair_type]][site_type] = {}
            adata[pair_type][cgenes[pair_type]][site_type]["present"] = None
            for cindex, _ in enumerate(comps.CLIP):
                adata[pair_type][cgenes[pair_type]][site_type]["clip%s" % cindex] = []

        stats[(proximal_reg, pair_type)] += 1
        stats[(proximal_reg, "all")] += 1
        gene_list.setdefault((proximal_reg, pair_type), []).append((gene_id, gene_name))
        gene_list.setdefault((proximal_reg, "all"), []).append((gene_id, gene_name))
        stats_bysite[("proximal", proximal_reg, pair_type)] += 1
        stats_bysite[("proximal", proximal_reg, "all")] += 1
        stats_bysite[("distal", distal_reg, pair_type)] += 1
        stats_bysite[("distal", distal_reg, "all")] += 1
        stats_bysite[("s1", s1_reg, pair_type)] += 1
        stats_bysite[("s1", s1_reg, "all")] += 1
        stats_bysite[("s2", s2_reg, pair_type)] += 1
        stats_bysite[("s2", s2_reg, "all")] += 1

        (_, proximal_lenup, proximal_lendown), (_, distal_lenup, distal_lendown), (_, s1_lenup, s1_lendown), (_, s2_lenup, s2_lendown) = coords(strand, proximal_pos, distal_pos, s1_pos, s2_pos)

        proximal_seq = pybio.genomes.seq(genome, chr, strand, proximal_pos, start=-proximal_lenup, stop=proximal_lendown)
        distal_seq = pybio.genomes.seq(genome, chr, strand, distal_pos, start=-distal_lenup, stop=distal_lendown)

        if s1_pos!=None:
            s1_seq = pybio.genomes.seq(genome, chr, strand, s1_pos, start=-s1_lenup, stop=s1_lendown)
        else:
            s1_seq = "N"
        if s2_pos!=None:
            s2_seq = pybio.genomes.seq(genome, chr, strand, s2_pos, start=-s2_lenup, stop=s2_lendown)
        else:
            s2_seq = "N"

        # count presence
        proximal_pre = presence_vector(proximal_seq, proximal_lenup, proximal_lendown, surr)
        distal_pre = presence_vector(distal_seq, distal_lenup, distal_lendown, surr)
        s1_pre = presence_vector(s1_seq, s1_lenup, s1_lendown, surr)
        s2_pre = presence_vector(s2_seq, s2_lenup, s2_lendown, surr)

        proximal_seq = adjust_len(proximal_seq, proximal_lenup, proximal_lendown, surr)
        distal_seq = adjust_len(distal_seq, distal_lenup, distal_lendown, surr)
        s1_seq = adjust_len(s1_seq, s1_lenup, s1_lendown, surr)
        s2_seq = adjust_len(s2_seq, s2_lenup, s2_lendown, surr)

        adata[pair_type][cgenes[pair_type]]["proximal"]["present"] = proximal_pre; adata[pair_type][cgenes[pair_type]]["distal"]["present"] = distal_pre; adata[pair_type][cgenes[pair_type]]["s1"]["present"] = s1_pre; adata[pair_type][cgenes[pair_type]]["s2"]["present"] = s2_pre

        present[("proximal", proximal_reg, pair_type)] = [x+y for x,y in zip(present[("proximal", proximal_reg, pair_type)], proximal_pre)]
        present[("proximal", proximal_reg, "all")] = [x+y for x,y in zip(present[("proximal", proximal_reg, "all")], proximal_pre)]
        present[("distal", distal_reg, pair_type)] = [x+y for x,y in zip(present[("distal", distal_reg, pair_type)], distal_pre)]
        present[("distal", distal_reg, "all")] = [x+y for x,y in zip(present[("distal", distal_reg, "all")], distal_pre)]
        present[("s1", s1_reg, pair_type)] = [x+y for x,y in zip(present[("s1", s1_reg, pair_type)], s1_pre)]
        present[("s1", s1_reg, "all")] = [x+y for x,y in zip(present[("s1", s1_reg, "all")], s1_pre)]
        present[("s2", s2_reg, pair_type)] = [x+y for x,y in zip(present[("s2", s2_reg, pair_type)], s2_pre)]
        present[("s2", s2_reg, "all")] = [x+y for x,y in zip(present[("s2", s2_reg, "all")], s2_pre)]

        if proximal_reg in ["enhanced", "repressed"]:
            fasta_files[("proximal", pair_type, proximal_reg)].write(">%s:%s %s%s:%s\n%s\n" % (gene_id, gene_name, strand, chr, proximal_pos, proximal_seq))
            fasta_files[("proximal", "all", proximal_reg)].write(">%s:%s %s%s:%s\n%s\n" % (gene_id, gene_name, strand, chr, proximal_pos, proximal_seq))
            fasta_files[("distal", pair_type, distal_reg)].write(">%s:%s %s%s:%s\n%s\n" % (gene_id, gene_name, strand, chr, distal_pos, distal_seq))
            fasta_files[("distal", "all", distal_reg)].write(">%s:%s %s%s:%s\n%s\n" % (gene_id, gene_name, strand, chr, distal_pos, distal_seq))

        # CLIP
        for cindex, clip_name in enumerate(comps.CLIP):
            for (site, reg, pos, len_up, len_down) in [("proximal", proximal_reg, proximal_pos, proximal_lenup, proximal_lendown), ("distal", distal_reg, distal_pos, distal_lenup, distal_lendown), ("s1", s1_reg, s1_pos, s1_lenup, s1_lendown), ("s2", s2_reg, s2_pos, s2_lenup, s2_lendown)]:
                adata[pair_type][cgenes[pair_type]][site]["reg"] = reg
                z = []
                if pos!=None:
                    z = clip[clip_name].get_vector("chr"+chr, strand, pos, -len_up, len_down)
                else:
                    z = [0]
                z = adjust_len(z, len_up, len_down, surr)
                if strand=="-":
                    z.reverse()
                cdata[(clip_name, site, reg, pair_type)] = [x+y for x,y in zip(cdata[(clip_name, site, reg, pair_type)], z)]
                cdata[(clip_name, site, reg, "all")] = [x+y for x,y in zip(cdata[(clip_name, site, reg, "all")], z)]
                # [(gene_id, gene_name), 0, 0, 1, 1 ....]
                z_vector = [(gene_id, gene_name, sum(z))] + z
                cdata_vectors[(clip_name, site, reg, pair_type)].append(z_vector)
                cdata_vectors[(clip_name, site, reg, "all")].append(z_vector)
                assert(len(cdata[(clip_name, site, reg, pair_type)])==401)
                assert(len(cdata[(clip_name, site, reg, "all")])==401)
                adata[pair_type][cgenes[pair_type]][site]["clip%s" % cindex] = z

        cgenes[pair_type] += 1
        r = f.readline()
    f.close() # end of reading gene data

    # save all data

    #print cdata[("peaks_id80654_rnd100_flank3_fdr0.05_group_5207_TARDBP-LC-flag-GFP-IP-group_sum_S_hg19--ensembl59_from_5158-5159_bedGraph-cDNA.bed.gz_lowFDR_clusters.bed.gz", "proximal", "r", "tandem")][:30]
    #print present[("proximal", "r", "tandem")][:30]

    temp = {}
    f = open(os.path.join(apa.path.comps_folder, comps_id, "%s_pair_distances.tab" % comps_id), "wt")
    dtemp = [0, 100, 450, 1000, 5000, 10000]
    for dist, count in pair_dist_stat.items():
        for x1, x2 in zip(dtemp, dtemp[1:]):
            if x1<=dist<x2:
                k = "%s..%s" % (x1, x2)
                temp[k] = temp.get(k, 0) + count
        if dist>10000:
            temp[">10000"] = temp.get(">10000", 0) + count
    for leq in ["0..100", "100..450", "450..1000", "1000..5000", "5000..10000", ">10000"]:
        f.write("%s nt\t%s pairs\n" % (leq, temp.get(leq, 0)))
    f.close()

    print "save adata.pickle"
    f = open(os.path.join(rnamap_dest, "adata.pickle"), "wb")
    f.write(cPickle.dumps(adata, protocol=-1))
    f.close()

    f = open(os.path.join(rnamap_dest, "stats.pickle"), "wb")
    f.write(cPickle.dumps(stats, protocol=-1))
    f.close()

    f = open(os.path.join(rnamap_dest, "stats_bysite.pickle"), "wb")
    f.write(cPickle.dumps(stats_bysite, protocol=-1))
    f.close()

    present_pairs = list(set([pair_type for (reg, pair_type) in stats.keys()]))
    sorted_present_pairs = []
    if "tandem" in present_pairs:
        sorted_present_pairs.append("tandem")
    if "composite" in present_pairs:
        sorted_present_pairs.append("composite")
    if "skipped" in present_pairs:
        sorted_present_pairs.append("skipped")
    if "all" in present_pairs:
        sorted_present_pairs.append("all")
    present_pairs = sorted_present_pairs

    #print stats.items()
    print "---"
    #print stats_bysite.items()
    for cat in [("proximal", "repressed", "tandem"), ("proximal", "enhanced", "tandem"), ("proximal", "repressed", "composite"), ("proximal", "enhanced", "composite"), ("proximal", "repressed", "skipped"), ("proximal", "enhanced", "skipped")]:
        print cat, stats_bysite.get(cat, 0)
    print "---"

    for f in fasta_files.values():
        f.close()

    # normalize with nt resolution
    print "normalization with nucleotide resolution"
    for pair_type in present_pairs:
        for site in ["proximal", "distal", "s1", "s2"]:
            for reg in ["enhanced", "repressed", "control_up", "control_down"]:
                n = present[(site, reg, pair_type)]
                if reg not in ["control_up", "control_down"]:
                    sdata[(site, reg, pair_type)] = [e/max(1.0, float(z)) for e,z in zip(sdata[(site, reg, pair_type)], n)]
                    pdata[(site, reg, pair_type)] = [e/max(1.0, float(z)) for e,z in zip(pdata[(site, reg, pair_type)], n)]
                for clip_name in comps.CLIP:
                    cdata[(clip_name, site, reg, pair_type)] = [e/max(1.0, float(z)) for e,z in zip(cdata[(clip_name, site, reg, pair_type)], n)]

    print "getting max values"
    cmax = {}
    fmax = {}
    for clip_name in comps.CLIP:
        cmax[clip_name] = {"tandem":0, "composite":0, "skipped":0, "all":0}
        fmax[clip_name] = {"tandem":0, "composite":0, "skipped":0, "all":0}
    smax = {"tandem":0, "composite":0, "skipped":0, "all":0}
    pmax = {"tandem":0, "composite":0, "skipped":0, "all":0}
    for pair_type in present_pairs:
        for site in ["proximal", "distal", "s1", "s2"]:
            for reg in ["enhanced", "repressed", "control_up", "control_down"]:
                for clip_name in comps.CLIP:
                    cmax[clip_name][pair_type] = max(cmax[clip_name][pair_type], max(cdata[(clip_name, site, reg, pair_type)]))
                    fmax[clip_name][pair_type] = max(fmax[clip_name][pair_type], rnamap_freq(cdata_vectors[(clip_name, site, "enhanced", pair_type)], cdata_vectors[(clip_name, site, "repressed", pair_type)], cdata_vectors[(clip_name, site, "control_up", pair_type)], cdata_vectors[(clip_name, site, "control_down", pair_type)], return_ymax=True))
                smax[pair_type] = max(smax[pair_type], max(sdata[(site, reg, pair_type)]))
                pmax[pair_type] = max(pmax[pair_type], max(pdata[(site, reg, pair_type)]))

    print "saving figures"
    for pair_type in present_pairs:
        for site in ["proximal", "distal", "s1", "s2"]:
            # clip
            if len(comps.CLIP)>0:
                for clip_index, clip_name in enumerate(comps.CLIP):
                    rnamap_area(cdata[(clip_name, site, "enhanced", pair_type)], cdata[(clip_name, site, "repressed", pair_type)], cdata[(clip_name, site, "control_up", pair_type)], cdata[(clip_name, site, "control_down", pair_type)], os.path.join(rnamap_dest, "clip%s.%s.%s" % (clip_index, pair_type, site)), title="%s.%s" % (pair_type, site), ymax=cmax[clip_name][pair_type], site=site, pair_type=pair_type, stats=stats)
                    rnamap_freq(cdata_vectors[(clip_name, site, "enhanced", pair_type)], cdata_vectors[(clip_name, site, "repressed", pair_type)], cdata_vectors[(clip_name, site, "control_up", pair_type)], cdata_vectors[(clip_name, site, "control_down", pair_type)], os.path.join(rnamap_dest, "clip%s_freq.%s.%s" % (clip_index, pair_type, site)), pair_type=pair_type, stats=stats, site=site, ymax=fmax[clip_name][pair_type])
                    rnamap_heat(cdata_vectors[(clip_name, site, "enhanced", pair_type)], cdata_vectors[(clip_name, site, "repressed", pair_type)], os.path.join(rnamap_dest, "clip%s_heat.%s.%s" % (clip_index, pair_type, site)), pair_type=pair_type, stats=stats, site=site, title=comps_id)
                    save_control_binding(cdata_vectors[(clip_name, site, "control_up", pair_type)], cdata_vectors[(clip_name, site, "control_down", pair_type)], os.path.join(rnamap_dest, "clip%s_heat.%s.%s" % (clip_index, pair_type, site)), site=site, pair_type=pair_type)

    # save gene lists (.tab files)
    for pair_type in present_pairs:
        f = open(os.path.join(rnamap_dest, "data_%s.tab" % pair_type), "wt")
        header = ["gene_id", "gene_name", "pair_type", "proximal_reg", "distal_reg", "pc", "fisher"]
        f.write("\t".join(header) + "\n")
        rev = {"enhanced":"repressed", "repressed":"enhanced", "control_up":"control_down", "control_down":"control_up", None:None}
        for proximal_reg in ["enhanced", "repressed", "control_up", "control_down"]:
            if gene_list.get((proximal_reg, pair_type), None)!=None:
                for gene_id, gene_name in gene_list[(proximal_reg, pair_type)]:
                    f.write("\t".join([gene_id, gene_name, pair_type, proximal_reg, rev[proximal_reg]])+"\n")
        # write python list of gene names
        if gene_list.get(("repressed", pair_type), None)!=None:
            L = ['"%s"' % gene_name for gene_id, gene_name in gene_list[("repressed", pair_type)]]
            f.write("genes_repressed=[%s]" % (",".join(L))+"\n")
        if gene_list.get(("enhanced", pair_type), None)!=None:
            L = ['"%s"' % gene_name for gene_id, gene_name in gene_list[("enhanced", pair_type)]]
            f.write("genes_enhanced=[%s]" % (",".join(L))+"\n")
        if gene_list.get(("control_up", pair_type), None)!=None:
            L = ['"%s"' % gene_name for gene_id, gene_name in gene_list[("control_up", pair_type)]]
            f.write("genes_controls_up=[%s]" % (",".join(L))+"\n")
        if gene_list.get(("control_down", pair_type), None)!=None:
            L = ['"%s"' % gene_name for gene_id, gene_name in gene_list[("control_down", pair_type)]]
            f.write("genes_controls_down=[%s]" % (",".join(L))+"\n")
        f.close()

    f = open(os.path.join(rnamap_dest, "index.html"), "wt")
    f.write("<html>\n")

    head = """

<head>
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
</head>
"""

    head += "<title>" + comps_id + "</title>" + """

<script type="text/javascript" src="https://apa-db.org/software/highslide/highslide/highslide.js"></script>
<link rel="stylesheet" type="text/css" href="https://apa-db.org/software/highslide/highslide/highslide.css" />

<link href='https://fonts.googleapis.com/css?family=Open+Sans:300italic,400italic,600italic,700italic,800italic,400,300,600,700,800' rel='stylesheet' type='text/css'>

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

* {
  font-family: 'Open Sans', sans-serif;
  font-size: 11px;
}

.highslide img {
   border: 0px;
   outline: none;
}

a {
    text-decoration: none;
    color: blue;
}

a:visited {
    color: blue;
}

</style>
"""

    f.write(head+"\n")
    f.write("<body>\n")

    body = comps_id + "<br>"
    body += """
    <a href="#" class="show_hide">+ details</a> &nbsp;&nbsp; |  &nbsp;&nbsp;
    <a href="https://docs.google.com/drawings/d/1_m4iZ1c9YwKI-NOWMSCSdg2IzEGedj-NaMTIlHqirc0/edit?usp=sharing" target="_dgs">+ definition of genomic sites</a>

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
    pair distance at least = """ + str(comps.pair_dist) + """
    <br>
    """

    if len(comps.CLIP)>0:
        for index, clip_name in enumerate(comps.CLIP):
            body += "CLIP = %s" % os.path.basename(clip_name)
            if index==0:
                body += " (main file)"
            body += "<br>"

    body += """
    </div>
    <br>
    Control
    <div style="padding-left: 10px;">
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
    <div style="padding-left: 10px;">
    """

    body += exp_print_out(comps.test)

    body += """

    <br>
    Gene lists:
    </br>
    <div style="padding-left: 10px;">
    """
    for pair_type in present_pairs:
        fname = "data_%s.tab" % pair_type
        body += "<a href=" + fname + " target=_" + pair_type + ">" + pair_type + " gene list</a><br>"

    body += """
    </div>

    </div>
    </div>

    <br><br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>CLIP</b> dataset:&nbsp;
    <select id="cmb_clip" name="cmb_clip" style="width:300px;" onchange="javascript:change_clip();">
    """
    nocache = str(random.randint(0, 1000000))
    for index, clip_name in enumerate(comps.CLIP):
        body += "<option value=" + str(index) + ">" + clip_name
        if index==0:
            body += " (main file)"
        body += "</option>"
    body += """

    </select><br><br>

<script>
function change_clip()
{
  index = $('#cmb_clip option:selected').val();
  sites = ["tandem", "composite", "skipped", "all"];
  for (var i=0; i<sites.length; i++)
  {
"""

    body += '$("#" + sites[i]+"_c00_img").attr("src", "clip" + index + "." + sites[i] + ".proximal.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c00_link").attr("href", "clip" + index + "." + sites[i] + ".proximal.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c01_img").attr("src", "clip" + index + "." + sites[i] + ".distal.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c01_link").attr("href", "clip" + index + "." + sites[i] + ".distal.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c02_img").attr("src", "clip" + index + "." + sites[i] + ".s1.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c02_link").attr("href", "clip" + index + "." + sites[i] + ".s1.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c03_img").attr("src", "clip" + index + "." + sites[i] + ".s2.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c03_link").attr("href", "clip" + index + "." + sites[i] + ".s2.png?nocache=%s");' % nocache

    body += '$("#" + sites[i]+"_c10_img").attr("src", "clip" + index + "_freq." + sites[i] + ".proximal.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c10_link").attr("href", "clip" + index + "_freq." + sites[i] + ".proximal.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c11_img").attr("src", "clip" + index + "_freq." + sites[i] + ".distal.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c11_link").attr("href", "clip" + index + "_freq." + sites[i] + ".distal.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c12_img").attr("src", "clip" + index + "_freq." + sites[i] + ".s1.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c12_link").attr("href", "clip" + index + "_freq." + sites[i] + ".s1.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c13_img").attr("src", "clip" + index + "_freq." + sites[i] + ".s2.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c13_link").attr("href", "clip" + index + "_freq." + sites[i] + ".s2.png?nocache=%s");' % nocache

    body += '$("#" + sites[i]+"_c20_img").attr("src", "clip" + index + "_heat." + sites[i] + ".proximal_pos.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c20_link").attr("href", "clip" + index + "_heat." + sites[i] + ".proximal_pos.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c21_img").attr("src", "clip" + index + "_heat." + sites[i] + ".distal_pos.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c21_link").attr("href", "clip" + index + "_heat." + sites[i] + ".distal_pos.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c22_img").attr("src", "clip" + index + "_heat." + sites[i] + ".s1_pos.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c22_link").attr("href", "clip" + index + "_heat." + sites[i] + ".s1_pos.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c23_img").attr("src", "clip" + index + "_heat." + sites[i] + ".s2_pos.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c23_link").attr("href", "clip" + index + "_heat." + sites[i] + ".s2_pos.png?nocache=%s");' % nocache

    body += '$("#" + sites[i]+"_c30_img").attr("src", "clip" + index + "_heat." + sites[i] + ".proximal_neg.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c30_link").attr("href", "clip" + index + "_heat." + sites[i] + ".proximal_neg.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c31_img").attr("src", "clip" + index + "_heat." + sites[i] + ".distal_neg.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c31_link").attr("href", "clip" + index + "_heat." + sites[i] + ".distal_neg.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c32_img").attr("src", "clip" + index + "_heat." + sites[i] + ".s1_neg.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c32_link").attr("href", "clip" + index + "_heat." + sites[i] + ".s1_neg.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c33_img").attr("src", "clip" + index + "_heat." + sites[i] + ".s2_neg.png?nocache=%s");' % nocache
    body += '$("#" + sites[i]+"_c33_link").attr("href", "clip" + index + "_heat." + sites[i] + ".s2_neg.png?nocache=%s");' % nocache

    body += """

  }
}
</script>

    """
    f.write(body+"\n")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px;'>")

    if len(comps.CLIP)>0:
        class_str = {"skipped":"skipped-exon", "composite":"composite-exon", "tandem":"same-exon", "all":"combined"}
        for t in present_pairs:
            ts = class_str[t]
            f.write("<tr><td align=center></td><td align=center>%s: proximal</td><td align=center>%s: distal</td><td align=center>%s: s1</td><td align=center>%s: s2</td></tr>\n" % (ts, ts, ts, ts))
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>iCLIP<br>enh=%s, rep=%s, con=%s</td>" % (ts, stats[("enhanced", t)], stats[("repressed", t)], stats[("control_up", t)]+stats[("control_down", t)]))
            f.write("<td align=right valign=center><a id=%s_c00_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c00_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0.%s.proximal.png?nocache=%s" % (t, nocache), t, "clip0.%s.proximal.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c01_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c01_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0.%s.distal.png?nocache=%s" % (t, nocache), t, "clip0.%s.distal.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c02_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c02_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0.%s.s1.png?nocache=%s" % (t, nocache), t, "clip0.%s.s1.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c03_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c03_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0.%s.s2.png?nocache=%s" % (t, nocache), t, "clip0.%s.s2.png?nocache=%s" % (t, nocache)))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>iCLIP<br>targets</td>" % ts)
            f.write("<td align=right valign=center><a id=%s_c10_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c10_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_freq.%s.proximal.png?nocache=%s" % (t, nocache), t, "clip0_freq.%s.proximal.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c11_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c11_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_freq.%s.distal.png?nocache=%s" % (t, nocache), t, "clip0_freq.%s.distal.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c12_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c12_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_freq.%s.s1.png?nocache=%s" % (t, nocache), t, "clip0_freq.%s.s1.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c13_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c13_img src=%s onerror=\"this.style.display='none'\"width=300px></a></td>" % (t, "clip0_freq.%s.s2.png?nocache=%s" % (t, nocache), t, "clip0_freq.%s.s2.png?nocache=%s" % (t, nocache)))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>iCLIP<br>enh</td>" % ts)
            f.write("<td align=right valign=center><a id=%s_c20_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c20_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_heat.%s.proximal_pos.png?nocache=%s" % (t, nocache), t, "clip0_heat.%s.proximal_pos.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c21_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c21_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_heat.%s.distal_pos.png?nocache=%s" % (t, nocache), t, "clip0_heat.%s.distal_pos.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c22_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c22_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_heat.%s.s1_pos.png?nocache=%s" % (t, nocache), t, "clip0_heat.%s.s1_pos.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c23_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c23_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_heat.%s.s2_pos.png?nocache=%s" % (t, nocache), t, "clip0_heat.%s.s2_pos.png?nocache=%s" % (t, nocache)))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>iCLIP<br>rep</td>" % ts)
            f.write("<td align=right valign=center><a id=%s_c30_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c30_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_heat.%s.proximal_neg.png?nocache=%s" % (t, nocache), t, "clip0_heat.%s.proximal_neg.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c31_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c31_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_heat.%s.distal_neg.png?nocache=%s" % (t, nocache), t, "clip0_heat.%s.distal_neg.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c32_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c32_img src=%s onerror=\"this.style.display='none'\" width=300px></a></td>" % (t, "clip0_heat.%s.s1_neg.png?nocache=%s" % (t, nocache), t, "clip0_heat.%s.s1_neg.png?nocache=%s" % (t, nocache)))
            f.write("<td align=right valign=center><a id=%s_c33_link href=%s class='highslide' onclick='return hs.expand(this)'><img id=%s_c33_img src=%s onerror=\"this.style.display='none'\"width=300px></a></td>" % (t, "clip0_heat.%s.s2_neg.png?nocache=%s" % (t, nocache), t, "clip0_heat.%s.s2_neg.png?nocache=%s" % (t, nocache)))
            f.write("</tr>")
            f.write("<tr><td><br><br></td><td><br><br></td><td><br><br></td><td><br><br></td><td><br><br></td></tr>")
            f.write("\n")

    f.write("\n")
    f.write("</table>")

    f.write("</body>")
    f.write("</html>\n")
    f.close()
