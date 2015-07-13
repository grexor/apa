import matplotlib
matplotlib.use("Agg", warn=False)
import matplotlib.pyplot as plt
import math
import gzip
from matplotlib import cm as CM
import numpy
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

def area_apa(vpos, vneg, filename, title="test", ymax=None):
    import matplotlib
    matplotlib.use("Agg", warn=False)
    import matplotlib.pyplot as plt
    import math
    import gzip
    from matplotlib import cm as CM
    import numpy
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
    a = plt.axes([0.1, 0.2, 0.8, 0.7])
    a.set_xlim(0, 400)
    if ymax!=None:
        a.set_ylim(-ymax, ymax)
    plt.ylabel("cpm divided by #genes")
    plt.xlabel("distance from poly-A site (nt)")

    vpos = pybio.utils.smooth(vpos)
    vneg = [-el for el in pybio.utils.smooth(vneg)]

    plt.fill_between(range(0, len(vpos)), 0, vpos, facecolor='red', alpha=0.6, interpolate=True)
    plt.plot(range(0, len(vpos)), vpos, color='red', alpha=1)

    plt.fill_between(range(0, len(vneg)), 0, vneg, facecolor='blue', alpha=0.6, interpolate=True)
    plt.plot(range(0, len(vneg)), vneg, color='blue', alpha=1)

    p = mpatches.Rectangle([200, -100], 0.01, 200, facecolor='none', edgecolor=(0.8, 0, 0))
    plt.gca().add_patch(p)
    plt.xticks([0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400], [-200,-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200])

    plt.title(title)
    plt.savefig(filename, dpi=100)
    print "saving", filename
    plt.close()

def process(comps_id=None, tab_file=None, clip_file="", genome=None, rnamap_dest="."):

    if comps_id!=None:
        comps = apa.comps.read_comps(comps_id)
        genome = comps.species
        if comps.iCLIP_filename!=None:
            clip_file = os.path.join(apa.path.iCLIP_folder, comps.iCLIP_filename)
        tab_file = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
        rnamap_dest = os.path.join(apa.path.comps_folder, comps_id, "rnamap")

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
    cdata = {}
    sdata = {}
    for pair_type in ["tandem", "composite", "skipped"]:
        cdata["e.%s.siteup" % pair_type] = [0] * 401
        cdata["e.%s.sitedown" % pair_type] = [0] * 401
        cdata["r.%s.siteup" % pair_type] = [0] * 401
        cdata["r.%s.sitedown" % pair_type] = [0] * 401
        sdata["e.%s.siteup" % pair_type] = [0] * 401
        sdata["e.%s.sitedown" % pair_type] = [0] * 401
        sdata["r.%s.siteup" % pair_type] = [0] * 401
        sdata["r.%s.sitedown" % pair_type] = [0] * 401
    f = open(tab_file, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        chr = data["chr"]
        strand = data["strand"]
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

        # UG
        # siteup
        _, z = pybio.sequence.search(seq_up, "TGTG")
        z = pybio.sequence.filter(z, hw=25, hwt=10)
        sdata["%s.%s.siteup" % (reg_siteup, pair_type)] = [x+y for x,y in zip(sdata["%s.%s.siteup" % (reg_siteup, pair_type)],z)]
        assert(len(sdata["%s.%s.siteup" % (reg_siteup, pair_type)])==401)

        # sitedown
        _, z = pybio.sequence.search(seq_down, "TGTG")
        z = pybio.sequence.filter(z, hw=25, hwt=10)
        sdata["%s.%s.sitedown" % (reg_sitedown, pair_type)] = [x+y for x,y in zip(sdata["%s.%s.sitedown" % (reg_sitedown, pair_type)],z)]
        assert(len(sdata["%s.%s.sitedown" % (reg_sitedown, pair_type)])==401)

        # CLIP
        # siteup
        z = []
        for index, x in enumerate(range(siteup_pos-200, siteup_pos+201)):
            z.append(clip.get_value("chr"+chr, strand, x, db="cpm"))
        if strand=="-":
            z.reverse()
        cdata["%s.%s.siteup" % (reg_siteup, pair_type)] = [x+y for x,y in zip(cdata["%s.%s.siteup" % (reg_siteup, pair_type)],z)]
        assert(len(cdata["%s.%s.siteup" % (reg_siteup, pair_type)])==401)

        # sitedown
        z = []
        for index, x in enumerate(range(sitedown_pos-200, sitedown_pos+201)):
            z.append(clip.get_value("chr" + chr, strand, x, db="cpm"))
        if strand=="-":
            z.reverse()
        cdata["%s.%s.sitedown" % (reg_sitedown, pair_type)] = [x+y for x,y in zip(cdata["%s.%s.sitedown" % (reg_sitedown, pair_type)],z)]
        assert(len(cdata["%s.%s.sitedown" % (reg_sitedown, pair_type)])==401)

        r = f.readline()
    f.close()

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
        area_apa(sdata["e.%s.siteup" % pair_type], sdata["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "seq.%s.siteup.png" % pair_type), title="%s.siteup" % pair_type, ymax=smax[pair_type])
        area_apa(sdata["e.%s.sitedown" % pair_type], sdata["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "seq.%s.sitedown.png" % pair_type), title="%s.sitedown" % pair_type, ymax=smax[pair_type])
        area_apa(cdata["e.%s.siteup" % pair_type], cdata["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "clip.%s.siteup.png" % pair_type), title="%s.siteup" % pair_type, ymax=cmax[pair_type])
        area_apa(cdata["e.%s.sitedown" % pair_type], cdata["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "clip.%s.sitedown.png" % pair_type), title="%s.sitedown" % pair_type, ymax=cmax[pair_type])

        # also write svg files
        area_apa(sdata["e.%s.siteup" % pair_type], sdata["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "seq.%s.siteup.svg" % pair_type), title="%s.siteup" % pair_type, ymax=smax[pair_type])
        area_apa(sdata["e.%s.sitedown" % pair_type], sdata["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "seq.%s.sitedown.svg" % pair_type), title="%s.sitedown" % pair_type, ymax=smax[pair_type])
        area_apa(cdata["e.%s.siteup" % pair_type], cdata["r.%s.siteup" % pair_type], os.path.join(rnamap_dest, "clip.%s.siteup.svg" % pair_type), title="%s.siteup" % pair_type, ymax=cmax[pair_type])
        area_apa(cdata["e.%s.sitedown" % pair_type], cdata["r.%s.sitedown" % pair_type], os.path.join(rnamap_dest, "clip.%s.sitedown.svg" % pair_type), title="%s.sitedown" % pair_type, ymax=cmax[pair_type])

    f = open(os.path.join(rnamap_dest, "index.html"), "wt")
    f.write("<html>\n")

    head = """<head>
    </head>"""

    f.write(head+"\n")

    f.write("<body>\n")

    body = """<div style="font-size: 12px;">
    <b>Description</b>
    <div style="font-size: 12px; padding-left: 10px;">
    <font color=red>red = pc>0</font>
    <br>
    <font color=blue>blue = pc<0</font>
    <br>
    pc = [ control(up_site) / sum(control) ] - [ test(up_site) / sum(test) ]
    </div>
    <br>
    <b>Parameters</b>
    <div style="font-size: 12px; padding-left: 10px;">
    pc threshold = """ + str(pc_thr) + """
    <br>
    fisher threshold = """ + str(fisher_thr) + """
    <br>
    pair distance at least = """ + str(pair_dist) + """
    <br>
    </div>
    </div>
    <br>
    <br>
    """

    f.write(body+"\n")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 12px;'>")
    f.write("<tr><td colspan=3 align=center style='font-size:13px'><b>iCLIP = " + os.path.basename(clip_file) + "</b><br><br></td></tr>")
    f.write("<tr><td align=center> </td><td align=center><b>siteup</td><td align=center><b>sitedown</td></tr>\n")
    f.write("<tr>")
    f.write("<td align=center valign=center><b>tandem</b><br><font color=red>e=%s</font><br><font color=blue>r=%s</font></b></td>" % (stats["e.tandem"], stats["r.tandem"]))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("clip.tandem.siteup.png", "clip.tandem.siteup.png"))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("clip.tandem.sitedown.png", "clip.tandem.sitedown.png"))
    f.write("</tr>")
    f.write("\n")

    f.write("<tr>")
    f.write("<td align=center valign=center><b>composite</b><br><font color=red>e=%s</font><br><font color=blue>r=%s</font></b></td>" % (stats["e.composite"], stats["r.composite"]))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("clip.composite.siteup.png", "clip.composite.siteup.png"))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("clip.composite.sitedown.png", "clip.composite.sitedown.png"))
    f.write("</tr>")
    f.write("\n")

    f.write("<tr>")
    f.write("<td align=center valign=center><b>skipped</b><br><font color=red>e=%s</font><br><font color=blue>r=%s</font></b></td>" % (stats["e.skipped"], stats["r.skipped"]))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("clip.skipped.siteup.png", "clip.skipped.siteup.png"))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("clip.skipped.sitedown.png", "clip.skipped.sitedown.png"))
    f.write("</tr>")
    f.write("\n")

    f.write("\n")
    f.write("</table>")
    f.write("<br><br><br>")

    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 12px;'>")
    f.write("<tr><td colspan=3 align=center><b>UG sequence search<br><br></b></td></tr>")
    f.write("<tr><td align=center> </td><td align=center><b>siteup</td><td align=center><b>sitedown</td></tr>\n")
    f.write("<tr>")
    f.write("<td align=center valign=center><b>tandem</b><br><font color=red>e=%s</font><br><font color=blue>r=%s</font></b></td>" % (stats["e.tandem"], stats["r.tandem"]))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("seq.tandem.siteup.png", "seq.tandem.siteup.png"))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("seq.tandem.sitedown.png", "seq.tandem.sitedown.png"))
    f.write("</tr>")
    f.write("\n")

    f.write("<tr>")
    f.write("<td align=center valign=center><b>composite</b><br><font color=red>e=%s</font><br><font color=blue>r=%s</font></b></td>" % (stats["e.composite"], stats["r.composite"]))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("seq.composite.siteup.png", "seq.composite.siteup.png"))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("seq.composite.sitedown.png", "seq.composite.sitedown.png"))
    f.write("</tr>")
    f.write("\n")

    f.write("<tr>")
    f.write("<td align=center valign=center><b>skipped</b><br><font color=red>e=%s</font><br><font color=blue>r=%s</font></b></td>" % (stats["e.skipped"], stats["r.skipped"]))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("seq.skipped.siteup.png", "seq.skipped.siteup.png"))
    f.write("<td align=right valign=center><a href=%s target=_new><img src=%s width=500px></a></td>" % ("seq.skipped.sitedown.png", "seq.skipped.sitedown.png"))
    f.write("</tr>")
    f.write("\n")

    f.write("\n")
    f.write("</table>")

    f.write("</body>")
    f.write("</html>\n")
    f.close()
