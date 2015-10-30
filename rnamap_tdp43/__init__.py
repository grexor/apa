# https://docs.google.com/document/d/1ZDxnRLuQWThZnd9CjfNSTyPmoIccrp4icaNXYkk6Y4s/edit#

import apa
import pybio
import numpy
import os
import numpy as np
from pandas import DataFrame
import math

def rnamap_area(vpos, vneg, filename, title="test", ymax=None, site="proximal"):
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
    a.set_xlim(0, 200)
    if ymax!=None:
        a.set_ylim(-ymax, ymax)
    plt.ylabel("average C")
    plt.xlabel("distance from PAS (nt)")

    if site=="proximal":
        vpos_draw = pybio.utils.smooth(vpos)
        vneg_draw = [-el for el in pybio.utils.smooth(vneg)]
        plt.fill_between(range(0, len(vpos_draw)), 0, vpos_draw, facecolor='red', alpha=0.6, interpolate=True)
        plt.plot(range(0, len(vpos_draw)), vpos_draw, color='red', alpha=1)
        plt.fill_between(range(0, len(vneg_draw)), 0, vneg_draw, facecolor='blue', alpha=0.6, interpolate=True)
        plt.plot(range(0, len(vneg_draw)), vneg_draw, color='blue', alpha=1)

    # turn the graph around for distal sites
    if site=="distal":
        vpos_draw = pybio.utils.smooth(vneg)
        vneg_draw = [-el for el in pybio.utils.smooth(vpos)]
        plt.fill_between(range(0, len(vpos_draw)), 0, vpos_draw, facecolor='blue', alpha=0.6, interpolate=True)
        plt.plot(range(0, len(vpos_draw)), vpos_draw, color='blue', alpha=1)
        plt.fill_between(range(0, len(vneg_draw)), 0, vneg_draw, facecolor='red', alpha=0.6, interpolate=True)
        plt.plot(range(0, len(vneg_draw)), vneg_draw, color='red', alpha=1)

    p = mpatches.Rectangle([100, -100], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
    plt.gca().add_patch(p)
    plt.xticks([0,25,50,75,100,125,150,175,200], [-100,-75,-50,-25,0,25,50,75,100])

    print "saving", filename
    plt.title(title)
    plt.savefig(filename+".png", dpi=100)
    plt.close()

def rnamap_r(vpos, vneg, filename, title="test", ymax=None, site="proximal"):
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
    a.set_xlim(0, 200)
    if ymax!=None:
        a.set_ylim(-ymax, ymax)
    plt.ylabel("average R")
    plt.xlabel("distance from PAS (nt)")

    if site=="proximal":
        vpos_draw = pybio.utils.smooth(vpos)
        vneg_draw = [-el for el in pybio.utils.smooth(vneg)]
        plt.plot(range(0, len(vpos_draw)), vpos_draw, color='red', alpha=1)
        plt.plot(range(0, len(vneg_draw)), vneg_draw, color='blue', alpha=1)

    # turn the graph around for distal sites
    if site=="distal":
        vpos_draw = pybio.utils.smooth(vneg)
        vneg_draw = [-el for el in pybio.utils.smooth(vpos)]
        plt.plot(range(0, len(vpos_draw)), vpos_draw, color='blue', alpha=1)
        plt.plot(range(0, len(vneg_draw)), vneg_draw, color='red', alpha=1)

    p = mpatches.Rectangle([100, -100], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
    p = mpatches.Rectangle([0, 0], 400, ymax*0.001, facecolor='none', edgecolor=(0.8, 0.8, 0.8), linestyle='dotted')
    plt.gca().add_patch(p)
    plt.xticks([0,25,50,75,100,125,150,175,200], [-100,-75,-50,-25,0,25,50,75,100])

    print "saving", filename
    plt.title(title)
    plt.savefig(filename+".png", dpi=100)
    plt.close()

def rnamap_heat(vpos, vneg, filename, title="test", site="proximal", stats=None, pair_type="tandem", alpha=0.8):
    """
    Draw RNA heatmaps
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

        dn = d.ix[indices, 1:] # get only subset of rows
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
        ax.set_xticklabels([-100,-50,0,50,100], minor=False)

        #ax.grid(False)
        ax.set_xlim(0, 201)
        ax.set_ylim(0, len(indices))
        ax.invert_yaxis()
        ax = plt.gca()
        for t in ax.xaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        for t in ax.yaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        p = mpatches.Rectangle([100, -100], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
        plt.gca().add_patch(p)
        plt.title("%s (top 20)" % title)
        plt.subplots_adjust(left=0.07, right=0.99, top=0.90, bottom=0.05)
        cbar = fig.colorbar(heatmap, fraction=0.01, pad=0.01)
        print "saving %s" % (filename+"_%s.png" % reg_type)
        plt.savefig(filename+"_%s.png" % reg_type, dpi=150)
    return

set1 = ["ATGTA", "TACGT", "TGTAT", "TGCAT", "GTGTT", "AGTGT", "TATGT", "TTGTG", "GCATG", "TGTCT", "TGTGA", "GAGTG", "CTGTG", "TGAGT", "GCGTG", "GGTGT", "TGTGC", "ACGTG", "GTGTC", "ATGCG", "TGCGC", "ATGTG", "GTACG", "GTGTA", "GTATG", "TGCGT", "CGTGT", "GTGCG", "CGCGC", "CGCGT", "GCGCG", "TGTGT", "GTGTG"]
set2 = ["TGTAT", "GTGTT", "TATGT", "GCGTG", "TGTGC", "ATGTG", "GTACG", "GTGTA", "GTATG", "TGCGT", "CGTGT", "GTGCG", "CGCGT", "TGTGT", "GTGTG", "TGAAT", "GAATG"]
set3 = ["TGTAT", "TATGT", "ATGTG", "GTGTA", "GTATG", "TGCGT", "CGTGT", "GTGCG", "TGTGT", "GTGTG"]

set4 = ["ATGTA", "TTTTT", "TTTGT", "TGTTT", "GTTTG"]
set5 = ["TGCGT", "CGTGT", "GTGCG", "TGTGT", "GTGTG", "CGCGT", "CGCGC", "GCGCG"]

def compute_C(fname, x):
    f = pybio.data.Fasta(fname)
    Cacc = [0] * 200
    Racc = [0] * 200
    Cvec = []
    ngenes = 0
    Rvec = []
    while f.read():
        seq = f.sequence[100:300]
        ngenes += 1
        _, s1 = pybio.sequence.search(seq, set1)
        _, s3 = pybio.sequence.search(seq, set3)
        _, s4 = pybio.sequence.search(seq, set4)
        _, s5 = pybio.sequence.search(seq, set5)
        N1 = []
        N3 = []
        for index in range(0, len(s1)):
            if s1[index]==1:
                val = sum(s1[max(0, index-15):index+15+1])
                N1.append(val)
            else:
                N1.append(0)
            if s3[index]==1:
                val = sum(s3[max(0, index-15):index+15+1])
                N3.append(val)
            else:
                N3.append(0)

        # compute S
        max_distal = max(N3[100:150])
        if max_distal==0:
            Z1 = 0.7
        else:
            Z1 = max_distal / 5.0
        max_proximal = max(N3[50:100])
        if max_proximal==0:
            Z2 = 0.7
        else:
            Z2 = max_proximal / 5.0

        # apply max_distal and max_proximal
        C = [0]*200
        for i in range(0, 100):
            C[i] = min(1, N1[i]*N3[i]*Z1*x)
        for i in range(100, 200):
            C[i] = min(1, N1[i]*N3[i]*Z2*x)

        # compute R values
        R = []
        for index, c in enumerate(C):
            if c>0:
                M4 = sum(s4[max(0, index-15):index+15+1])
                M5 = sum(s5[max(0, index-15):index+15+1])
                if M4==0 and M5==0:
                    R.append(0)
                else:
                    Rval = c * math.log(float(1+M4)/(1+M5), 2)
                    R.append(Rval)
            else:
                R.append(0)

        assert(len(C)==len(R))

        Racc = [x1+y1 for x1,y1 in zip(R, Racc)]
        Cacc = [x1+y1 for x1,y1 in zip(C, Cacc)]
        Cvec.append([(f.id.split(" ")[0].split(":")[0], f.id.split(" ")[0].split(":")[1], sum(C))] + C)

    Cacc = [e/float(ngenes) for e in Cacc] # average
    Racc = [e/float(ngenes) for e in Racc] # average
    return ngenes, Cacc, Cvec, Racc

def process(comps_id):

    rnamap_folder = os.path.join(apa.path.comps_folder, comps_id, "rnamap_tdp43")
    if not os.path.exists(rnamap_folder):
        os.makedirs(rnamap_folder)

    #fasta_folder = os.path.join(apa.path.comps_folder, comps_id, "rnamap_pas")
    fasta_folder = os.path.join(apa.path.comps_folder, comps_id, "rnamap")

    xlist = [0.1, 0.01, 0.005, 0.003, 0.002, 0.001]

    stats = {}
    ymax = {}
    rymax = {}

    for x in xlist:
        for pair in ["tandem", "composite", "skipped"]:
            ymax["%s_%s" % (x, pair)] = 0
            rymax["%s_%s" % (x, pair)] = 0
            for site in ["proximal", "distal"]:
                stats["%s_%s_r" % (site, pair)], vneg, vnegv, rneg = compute_C(os.path.join(fasta_folder, "%s_%s_r.fasta" % (site, pair)), x=x)
                stats["%s_%s_e" % (site, pair)], vpos, vposv, rpos = compute_C(os.path.join(fasta_folder, "%s_%s_e.fasta" % (site, pair)), x=x)
                vneg = pybio.utils.smooth(vneg)
                vpos = pybio.utils.smooth(vpos)
                ymax["%s_%s" % (x, pair)] = max(ymax["%s_%s" % (x, pair)], max(vneg), max(vpos))
                rneg = pybio.utils.smooth(rneg)
                rpos = pybio.utils.smooth(rpos)
                rymax["%s_%s" % (x, pair)] = max(rymax["%s_%s" % (x, pair)], abs(max(rneg, key=abs)), abs(max(rpos, key=abs)))

    for x in xlist:
        for site in ["proximal", "distal"]:
            for pair in ["tandem", "composite", "skipped"]:
                print x
                stats["%s_%s_r" % (site, pair)], vneg, vnegv, rneg = compute_C(os.path.join(fasta_folder, "%s_%s_r.fasta" % (site, pair)), x=x)
                stats["%s_%s_e" % (site, pair)], vpos, vposv, rpos = compute_C(os.path.join(fasta_folder, "%s_%s_e.fasta" % (site, pair)), x=x)
                rnamap_area(vpos, vneg, os.path.join(rnamap_folder, "%s_%s_x%s" % (site, pair, x)), title="%s, %s, x=%s" % (site, pair, x), ymax=ymax["%s_%s" % (x, pair)], site=site)
                rnamap_heat(vposv, vnegv, os.path.join(rnamap_folder, "%s_%s_x%s_heat" % (site, pair, x)), title="%s, %s, x=%s" % (site, pair, x), site=site, stats=None, pair_type="tandem", alpha=0.8)
                rnamap_r(rpos, rneg, os.path.join(rnamap_folder, "%s_%s_x%s_r" % (site, pair, x)), title="%s, %s, x=%s" % (site, pair, x), site=site, ymax=rymax["%s_%s" % (x, pair)])

    f = open(os.path.join(rnamap_folder, "index.html"), "wt")
    f.write("<html>\n")

    head = """

    <head>
    </head>

    <script type="text/javascript" src="https://apa-db.org/software/highslide/highslide/highslide.js"></script>
    <link rel="stylesheet" type="text/css" href="https://apa-db.org/software/highslide/highslide/highslide.css" />

    <script type="text/javascript">
        hs.graphicsDir = 'https://apa-db.org/software/highslide/highslide/graphics/';
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

    body = """
    """
    f.write(body+"\n")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 12px;'>")

    for t in ["tandem", "composite", "skipped"]:
        f.write("<tr><td align=center></td><td align=center>proximal</td><td align=center>distal</td></tr>\n")
        for x in xlist:
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>x=%s<br>e=%s<br>r=%s</td>" % (t, x, stats["proximal_%s_e" % t], stats["proximal_%s_r" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("proximal_%s_x%s.png" % (t, x), "proximal_%s_x%s.png" % (t, x)))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("distal_%s_x%s.png" % (t, x), "distal_%s_x%s.png" % (t, x)))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>x=%s<br>e=%s<br>r=%s</td>" % (t, x, stats["proximal_%s_e" % t], stats["proximal_%s_r" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("proximal_%s_x%s_r.png" % (t, x), "proximal_%s_x%s_r.png" % (t, x)))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("distal_%s_x%s_r.png" % (t, x), "distal_%s_x%s_r.png" % (t, x)))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>x=%s<br>e=%s<br>r=%s</td>" % (t, x, stats["proximal_%s_e" % t], stats["proximal_%s_r" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("proximal_%s_x%s_heat_pos.png" % (t, x), "proximal_%s_x%s_heat_pos.png" % (t, x)))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("distal_%s_x%s_heat_pos.png" % (t, x), "distal_%s_x%s_heat_pos.png" % (t, x)))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>x=%s<br>e=%s<br>r=%s</td>" % (t, x, stats["proximal_%s_e" % t], stats["proximal_%s_r" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("proximal_%s_x%s_heat_neg.png" % (t, x), "proximal_%s_x%s_heat_neg.png" % (t, x)))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>" % ("distal_%s_x%s_heat_neg.png" % (t, x), "distal_%s_x%s_heat_neg.png" % (t, x)))
            f.write("</tr>")
            f.write("\n")
    f.write("\n")
    f.write("</table>")

    f.write("</body>")
    f.write("</html>\n")
    f.close()

def compute_C_new(fname, x, motif_set):
    f = pybio.data.Fasta(fname)
    Cacc = [0] * 200
    Racc = [0] * 200
    Cvec = []
    ngenes = 0
    Rvec = []
    while f.read():
        seq = f.sequence[100:300]
        ngenes += 1
        _, s1 = pybio.sequence.search(seq, motif_set)
        N1 = []
        for index in range(0, len(s1)):
            if s1[index]==1:
                val = sum(s1[max(0, index-15):index+15+1])
                N1.append(val)
            else:
                N1.append(0)

        # apply max_distal and max_proximal
        C = [0]*200
        for i in range(0, len(C)):
            assert(N1[i]<=31)
            C[i] = min(1, N1[i]*x)

        Cacc = [x1+y1 for x1,y1 in zip(C, Cacc)]
        Cvec.append([(f.id.split(" ")[0].split(":")[0], f.id.split(" ")[0].split(":")[1], sum(C))] + C)

    Cacc = [e/float(ngenes) for e in Cacc] # average
    return ngenes, Cacc, Cvec

def process2(comps_id, folder_name, motif_set):

    rnamap_folder = os.path.join(apa.path.comps_folder, comps_id, folder_name)
    if not os.path.exists(rnamap_folder):
        os.makedirs(rnamap_folder)

    #fasta_folder = os.path.join(apa.path.comps_folder, comps_id, "rnamap_pas")
    fasta_folder = os.path.join(apa.path.comps_folder, comps_id, "rnamap")

    #xlist = [0.1, 0.01, 0.005, 0.003, 0.002, 0.001]
    xlist = [0.12, 0.1, 0.07, 0.04]

    stats = {}
    ymax = {}

    for x in xlist:
        for pair in ["tandem", "composite", "skipped"]:
            ymax["%s_%s" % (x, pair)] = 0
            for site in ["proximal", "distal"]:
                stats["%s_%s_r" % (site, pair)], vneg, vnegv = compute_C_new(os.path.join(fasta_folder, "%s_%s_r.fasta" % (site, pair)), x=x, motif_set=motif_set)
                stats["%s_%s_e" % (site, pair)], vpos, vposv = compute_C_new(os.path.join(fasta_folder, "%s_%s_e.fasta" % (site, pair)), x=x, motif_set=motif_set)
                vneg = pybio.utils.smooth(vneg)
                vpos = pybio.utils.smooth(vpos)
                ymax["%s_%s" % (x, pair)] = max(ymax["%s_%s" % (x, pair)], max(vneg), max(vpos))

    for x in xlist:
        for site in ["proximal", "distal"]:
            for pair in ["tandem", "composite", "skipped"]:
                print x
                stats["%s_%s_r" % (site, pair)], vneg, vnegv = compute_C_new(os.path.join(fasta_folder, "%s_%s_r.fasta" % (site, pair)), x=x, motif_set=motif_set)
                stats["%s_%s_e" % (site, pair)], vpos, vposv = compute_C_new(os.path.join(fasta_folder, "%s_%s_e.fasta" % (site, pair)), x=x, motif_set=motif_set)
                rnamap_area(vpos, vneg, os.path.join(rnamap_folder, "%s_%s_x%s" % (site, pair, x)), title="%s, %s, x=%s" % (site, pair, x), ymax=ymax["%s_%s" % (x, pair)], site=site)
                rnamap_heat(vposv, vnegv, os.path.join(rnamap_folder, "%s_%s_x%s_heat" % (site, pair, x)), title="%s, %s, x=%s" % (site, pair, x), site=site, stats=None, pair_type="tandem", alpha=0.8)

    f = open(os.path.join(rnamap_folder, "index.html"), "wt")
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
      //$("#tandem_div").hide();
      $("#composite_div").hide();
      $("#skipped_div").hide();

      $('#tandem_btn').click(function (e) {
        $("#tandem_div").slideToggle("fast");
        var val = $(this).text() == "- tandem" ? "+ tandem" : "- tandem";
        $(this).hide().text(val).fadeIn("fast");
        e.preventDefault();
      });

      $('#composite_btn').click(function (e) {
        $("#composite_div").slideToggle("fast");
        var val = $(this).text() == "- composite" ? "+ composite" : "- composite";
        $(this).hide().text(val).fadeIn("fast");
        e.preventDefault();
      });

      $('#skipped_btn').click(function (e) {
        $("#skipped_div").slideToggle("fast");
        var val = $(this).text() == "- skipped" ? "+ skipped" : "- skipped";
        $(this).hide().text(val).fadeIn("fast");
        e.preventDefault();
      });

    });

    </script>

    <style>

    .show_hide {
      font-size: 16px;
    }

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

    body = """
    """
    f.write(body+"\n")

    for t in ["tandem", "composite", "skipped"]:

        if t=="tandem":
            f.write("<a href='#' class='show_hide' id='%s_btn'>- %s</a>\n" % (t,t)) # tandem is already open
        else:
            f.write("<a href='#' class='show_hide' id='%s_btn'>+ %s</a>\n" % (t,t)) # composite and skipped are closed

        f.write("<div id='%s_div' style='font-size: 12px; padding-left: 10px;' class='slidingDiv'>\n" % t)
        f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 12px;'>\n")
        f.write("<tr><td align=center></td><td align=center>proximal</td><td align=center>distal</td></tr>\n")
        for x in xlist:
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>x=%s<br>e=%s<br>r=%s</td>" % (t, x, stats["proximal_%s_e" % t], stats["proximal_%s_r" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>\n" % ("proximal_%s_x%s.png" % (t, x), "proximal_%s_x%s.png" % (t, x)))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>\n" % ("distal_%s_x%s.png" % (t, x), "distal_%s_x%s.png" % (t, x)))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>x=%s<br>e=%s<br>r=%s</td>" % (t, x, stats["proximal_%s_e" % t], stats["proximal_%s_r" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>\n" % ("proximal_%s_x%s_heat_pos.png" % (t, x), "proximal_%s_x%s_heat_pos.png" % (t, x)))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>\n" % ("distal_%s_x%s_heat_pos.png" % (t, x), "distal_%s_x%s_heat_pos.png" % (t, x)))
            f.write("</tr>")
            f.write("<tr>")
            f.write("<td align=right valign=center>%s<br>x=%s<br>e=%s<br>r=%s</td>" % (t, x, stats["proximal_%s_e" % t], stats["proximal_%s_r" % t]))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>\n" % ("proximal_%s_x%s_heat_neg.png" % (t, x), "proximal_%s_x%s_heat_neg.png" % (t, x)))
            f.write("<td align=right valign=center><a href=%s class='highslide' onclick='return hs.expand(this)'><img src=%s width=500px></a></td>\n" % ("distal_%s_x%s_heat_neg.png" % (t, x), "distal_%s_x%s_heat_neg.png" % (t, x)))
            f.write("</tr>")
            f.write("\n")
        f.write("\n")
        f.write("</table>\n")
        f.write("</div>\n")
        f.write("<br><br>\n")

    f.write("</body>\n")
    f.write("</html>\n")
    f.close()
