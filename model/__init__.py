import apa
import os
import sys
import pybio
from sklearn import svm
from sklearn import linear_model
import numpy as np
from sklearn import cross_validation
from sklearn import preprocessing
from sklearn.cross_validation import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
import itertools
import sklearn
from sklearn.metrics import roc_curve, auc
import pylab as pl
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import shutil
from sklearn import tree
import pickle
from scipy.interpolate import interp1d
from sklearn.metrics import explained_variance_score
from sklearn.metrics import mean_absolute_error
from sklearn.decomposition import PCA
#from sklearn.lda import LDA
import math
from collections import Counter
from sklearn.ensemble import ExtraTreesClassifier
import random

matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['axes.titlesize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['ytick.labelsize'] = 10
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['legend.frameon'] = 'False'

ntmap = {"A":0, "C":1, "T":2, "G":3}
cross_ratio = 0.4
seq_from = -50
seq_to = 50
len_seq = seq_to - seq_from + 1
clip_bin = 10

def make_motifs(size):
    motifs = []
    z = itertools.product("ATCG", repeat=size)
    for el in z:
        motifs.append("".join(el))
    return motifs

types = [("siteup_e", ["siteup_e", "sitedown_r"], ["siteup_c", "sitedown_c"]), ("siteup_r", ["siteup_r", "sitedown_e"], ["siteup_c", "sitedown_c"])]
types = [("siteup_e", ["siteup_e"], ["siteup_c"]), ("siteup_r", ["siteup_r"], ["siteup_c"]), ("sitedown_e", ["sitedown_e"], ["sitedown_c"]), ("sitedown_r", ["sitedown_r"], ["sitedown_c"])]
types = [("siteup", ["siteup_e"], ["siteup_r"]), ("sitedown", ["sitedown_e"], ["sitedown_r"])]
models = ["randomf", "svm"]
ftypes = ["A", "B", "C"]
ftypes = ["A", "D"]
#ftypes = ["A", "B", "D", "E", "F"]
ftypes = ["A", "D", "E", "F"]

k3 = {}
for index, m in enumerate(make_motifs(3)):
    k3[m] = index

k4 = {}
for index, m in enumerate(make_motifs(4)):
    k4[m] = index

k5 = {}
for index, m in enumerate(make_motifs(5)):
    k5[m] = index

def pwm_score(matrix, seq):
    s = 1
    for index, nt in enumerate(seq):
        s *= matrix[index][nt]
    return s

def run(comps_id, debug=False):
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")

    # do not prepare data every time when debugging, time consuming
    if not debug:
        if os.path.exists(model_folder):
            shutil.rmtree(model_folder)
        try:
            os.makedirs(model_folder)
        except:
            pass
        make_fasta(comps_id)
        prepare(comps_id)

    #plot_lda(comps_id)
    plot_pca(comps_id)

    # check
    #plot_hist(comps_id)

    # classification
    if "randomf" in models:
        predict_randomf(comps_id)
    if "svm" in models:
        predict_svm(comps_id)
    #predict_tree(comps_id)
    plot_roc(comps_id)
    write_index(comps_id)

    # regression
    #predict_linreg(comps_id)

# par = [pc, fisher, site_distance, site_type]
def make_fasta(comps_id):
    print "%s:model:make fasta" % comps_id
    comps = apa.comps.read_comps(comps_id)
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")

    # fasta
    fasta_folder = os.path.join(apa.path.comps_folder, comps_id, "fasta")
    if os.path.exists(fasta_folder):
        shutil.rmtree(fasta_folder)
    os.makedirs(fasta_folder)

    # rnamotifs
    rnamotifs_folder = os.path.join(apa.path.comps_folder, comps_id, "rnamotifs")
    if os.path.exists(rnamotifs_folder):
        shutil.rmtree(rnamotifs_folder)
    os.makedirs(rnamotifs_folder)

    # rnamotifs init
    for site in ["siteup", "sitedown"]:
        rnamotifs_filename = os.path.join(apa.path.comps_folder, comps_id, "rnamotifs", "%s_%s.tab" % (comps_id, site))
        rnamotifs_file = open(rnamotifs_filename, "wt")
        rnamotifs_file.write("\t".join(["id", "chr", "strand", "pos", "event_class"])+"\n")

    voom = {}
    voom_fname = os.path.join(apa.path.comps_folder, comps_id, "%s.voom_exons.tab" % comps_id)
    f = open(voom_fname, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        gene_id = r[2]
        exon_id = r[1]
        fdr = float(r[-1])
        voom["%s:%s" % (gene_id, exon_id)] = fdr
        r = f.readline()
    f.close()

    # r = repressed, e = enhanced, c = control
    stats = Counter()
    ntdist = {}
    for sitetype in ["siteup_e", "siteup_r", "sitedown_e", "sitedown_r", "siteup_c", "sitedown_c"]:
        ntdist[sitetype] = {}
        for index in range(0, len_seq):
            ntdist[sitetype][index] = {0:0, 1:0, 2:0, 3:0}
    input_file = os.path.join(apa.path.comps_folder, comps_id, "%s.pairs_de.tab" % comps_id)
    f = open(input_file, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    row_id = 0
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        chr = data["chr"]
        strand = data["strand"]
        gene_id = data["gene_id"]
        siteup_pos = int(data["siteup_pos"])
        sitedown_pos = int(data["sitedown_pos"])
        pc = float(data["pc"])
        fisher = float(data["fisher"])
        pair_type = data["pair_type"]
        fdr_voom_up = voom.get("%s:%s" % (gene_id, siteup_pos), 1)
        fdr_voom_down = voom.get("%s:%s" % (gene_id, sitedown_pos), 1)

        if abs(siteup_pos-sitedown_pos)<comps.pair_dist:
            r = f.readline()
            continue

        if pair_type!="tandem":
            r = f.readline()
            continue

        """
        # filtering by Fisher
        reg_siteup = None
        if pc>0 and abs(pc)>comps.pc_thr and fisher<comps.fisher_thr:
            reg_siteup = "e"
        if pc<0 and abs(pc)>comps.pc_thr and fisher<comps.fisher_thr:
            reg_siteup = "r"
        if abs(pc)<comps.control_thr:
            reg_siteup = "c"
        if reg_siteup==None:
            r = f.readline()
            continue
        """

        # filtering by voom
        reg_siteup = None
        if pc>0 and abs(pc)>comps.pc_thr and fdr_voom_up<comps.fisher_thr and fdr_voom_down<comps.fisher_thr:
            reg_siteup = "e"
        if pc<0 and abs(pc)>comps.pc_thr and fdr_voom_up<comps.fisher_thr and fdr_voom_down<comps.fisher_thr:
            reg_siteup = "r"
        if abs(pc)<comps.control_thr:
            reg_siteup = "c"
        if reg_siteup==None:
            r = f.readline()
            continue

        reg_sitedown = {"e":"r", "r":"e", "c":"c"}[reg_siteup]
        row_id += 1
        for (site_reg, site_type, site_pos) in [(reg_siteup, "siteup", siteup_pos), (reg_sitedown, "sitedown", sitedown_pos)]:
            reg_key = "%s_%s" % (site_type, site_reg)
            stats[reg_key] += 1
            seq = pybio.genomes.seq(comps.species, chr, strand, site_pos, start=seq_from, stop=seq_to)
            for index, n in enumerate(seq):
                ntdist[reg_key][index][ntmap[n]] += 1
            fasta_filename = os.path.join(apa.path.comps_folder, comps_id, "fasta", "%s.fasta" % (reg_key))
            if not os.path.exists(fasta_filename):
                fasta = open(fasta_filename, "wt")
            else:
                fasta = open(fasta_filename, "at")
            fasta.write(">%s.%s loci=%s:%s:%s:%s\n%s\n" % (reg_key, stats[reg_key], chr, strand, site_pos+seq_from, site_pos+seq_to, seq))
            fasta.close()

            # rnamotifs
            rnamotifs_filename = os.path.join(apa.path.comps_folder, comps_id, "rnamotifs", "%s_%s.tab" % (comps_id, site_type))
            rnamotifs_file = open(rnamotifs_filename, "at")
            rnamotifs_file.write("\t".join(str(e) for e in [row_id, chr, strand, site_pos, site_reg])+"\n")
            rnamotifs_file.close()

        r = f.readline()

    # plot nt-dist log ratio
    plot_data = {}
    ymax = 0
    for sitetype in ["siteup_e", "siteup_r", "sitedown_e", "sitedown_r", "siteup_c", "sitedown_c"]:
        sitecontrol = "siteup_c" if sitetype.startswith("siteup") else "sitedown_c"
        for nuc in [0, 1, 2, 3]:
            yt = [float(ntdist[sitetype][index][nuc])/(ntdist[sitetype][index][0]+ntdist[sitetype][index][1]+ntdist[sitetype][index][2]+ntdist[sitetype][index][3]) for index in range(0, len_seq)]
            yc = [float(ntdist[sitecontrol][index][nuc])/(ntdist[sitecontrol][index][0]+ntdist[sitecontrol][index][1]+ntdist[sitecontrol][index][2]+ntdist[sitecontrol][index][3]) for index in range(0, len_seq)]
            yt = [x*100 for x in yt]
            yc = [x*100 for x in yc]
            #y = [math.log((float(x+1)/(y+1)), 2) for x,y in zip(yt, yc)]
            y = [x-y for x,y in zip(yt, yc)]
            y = pybio.sequence.convolve(y, 1) # hw=1, sum with a window of 3
            y = [e/3.0 for e in y] # average
            ymax = max(ymax, max(y), min(y))
            plot_data["%s:%s" % (sitetype, nuc)] = y
    for sitetype in ["siteup_e", "siteup_r", "sitedown_e", "sitedown_r", "siteup_c", "sitedown_c"]:
        pl.clf()
        a = pl.axes([0.1, 0.1, 0.85, 0.8])
        a.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        a.spines["top"].set_alpha(0.1)
        a.spines["bottom"].set_alpha(0.1)
        a.spines["left"].set_alpha(0.1)
        a.spines["right"].set_alpha(0.1)
        for nuc in [0, 1, 2, 3]:
            pl.plot(range(0, len_seq), pybio.utils.smooth(plot_data["%s:%s" % (sitetype, nuc)]), label={0:"A", 1:"C", 2:"T", 3:"G"}[nuc])
        pl.ylim([-abs(ymax), abs(ymax)])
        pl.xticks([0,25,50,75,100], [-50,-25,0,25,50])
        pl.xlabel('distance [nt]')
        pl.ylabel('freq')
        pl.title('%s nucleotide distribution' % sitetype)
        pl.legend(loc="upper right")
        pl.savefig(os.path.join(model_folder, "%s_ntdist.png" % sitetype))
        pl.savefig(os.path.join(model_folder, "%s_ntdist.pdf" % sitetype))

    print stats.items()

def prepare(comps_id):

    def add_features(seq_rec, features, pos_from=None):
        seq, seq_id = seq_rec
        seq_id = seq_id.split("loci=")[1]
        seq_id = seq_id.split(":")
        chr = seq_id[0]
        strand = seq_id[1]
        pos_from = int(seq_id[2])
        pos_to = int(seq_id[3])
        # single positional nucleotides
        if "CLIP" in ftypes:
            for index in range(0, len_seq, clip_bin):
                features.append(iclip.get_region(chr, strand, int(pos_from)+index, int(pos_from)+index+clip_bin))

        if "A" in ftypes:
            for index in range(0, len_seq):
                features.append(ntmap[seq[index]])

        if "B" in ftypes:
            # 3-mers
            for index in range(0, len_seq):
                features.append(k3.get(seq[index:index+3], 0))

        if "C" in ftypes:
            # 4-mers
            for index in range(0, len_seq):
                features.append(k4.get(seq[index:index+5], 0))

        if "D" in ftypes:
            # count (no position)
            for m in make_motifs(3):
                features.append(seq[0:len(seq)/2].count(m))
                features.append(seq[len(seq)/2:].count(m))

        if "E" in ftypes:
            # count (no position)
            for m in make_motifs(4):
                features.append(seq[0:len(seq)/2].count(m))
                features.append(seq[len(seq)/2:].count(m))

        if "F" in ftypes:
            # count (no position)
            for m in make_motifs(5):
                features.append(seq[0:len(seq)/2].count(m))
                features.append(seq[len(seq)/2:].count(m))

        return features

    def construct_features(class_val, keys, seqs):
        result = []
        temp = []
        for k in keys:
            temp.append(seqs[k])
        for seq_index in range(0, len(temp[0])):
            features = [class_val]
            for r in range(0, len(temp)):
                features = add_features(temp[r][seq_index], features)
            result.append(features)
        return result

    comps = apa.comps.read_comps(comps_id)
    if comps.iCLIP_filename not in [None, "None", ""] and "CLIP" in ftypes:
        iclip_filename = os.path.join(apa.path.iCLIP_folder, comps.iCLIP_filename)
        iclip = pybio.data.Bedgraph(iclip_filename, genome=comps.species, fast=True) # we specify genome since iCLIP data has UCSC chromosome names, this converts it to Ensembl
    else:
        iclip = None

    print "%s:model:preparing data" % comps_id

    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    fasta_folder = os.path.join(apa.path.comps_folder, comps_id, "fasta")
    if not os.path.exists(model_folder):
        os.makedirs(model_folder)

    # types: (name, class_0, class_1)
    #types = [("siteup_e", ["siteup_e"], ["siteup_c", "sitedown_c"]), ("siteup_r", ["siteup_r"], ["siteup_c", "sitedown_c"])]

    for (name, class_0, class_1) in types:
        # name features
        features = []
        #if iclip!=None:
        #    for index in range(0, len_seq):
        #        features.append("clip_%s" % index)
        for site in class_0:
            if "CLIP" in ftypes:
                for index in range(0, len_seq, clip_bin):
                    features.append("%s_clip%s" % (site, (index+seq_from)))

            if "A"in ftypes:
                for index in range(0, len_seq):
                    features.append("%s_n%s" % (site, (index+seq_from)))

            if "B" in ftypes:
                for index in range(seq_from, seq_to+1):
                    features.append("k3_p%s" % index)

            if "C" in ftypes:
                for index in range(seq_from, seq_to+1):
                    features.append("k4_p%s" % index)

            if "D"in ftypes:
                search_motifs = make_motifs(3)
                for motif in search_motifs:
                    features.append("%s_c%s_up" % (site, motif))
                    features.append("%s_c%s_down" % (site, motif))

            if "E"in ftypes:
                search_motifs = make_motifs(4)
                for motif in search_motifs:
                    features.append("%s_c%s_up" % (site, motif))
                    features.append("%s_c%s_down" % (site, motif))

            if "F"in ftypes:
                search_motifs = make_motifs(5)
                for motif in search_motifs:
                    features.append("%s_c%s_up" % (site, motif))
                    features.append("%s_c%s_down" % (site, motif))

        # read in all sequences
        seqs = {}
        for site_name in class_0+class_1:
            seqs[site_name] = []
            fasta_file = os.path.join(fasta_folder, "%s.fasta" % site_name)
            f = pybio.data.Fasta(fasta_file)
            while f.read():
                seqs[site_name].append((f.sequence, f.id))

        data = []
        data.extend(construct_features(0, class_0, seqs))
        data.extend(construct_features(1, class_1, seqs))

        data_filename = os.path.join(model_folder, "%s.tab" % name)
        print data_filename
        header = ["id", "class"]
        header.extend(features)
        fr = open(data_filename, "wt")
        fr.write("\t".join(str(x) for x in header) + "\n")
        id = 1
        for rec in data:
            rec.insert(0, id)
            id += 1
            fr.write("\t".join(str(x) for x in rec) + "\n")
        fr.close()

def predict_svm(comps_id):
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    for fn, _, _ in types:
        print "%s:model:svm.%s" % (comps_id, fn)
        data_file = open(os.path.join(model_folder, "%s.tab" % fn), "rt")
        header = data_file.readline().replace("\r", "").replace("\n", "").split("\t")
        r = data_file.readline()
        y = []
        x = []
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            x_row = []
            for el in r[3:]:
                x_row.append(float(el))
            x.append(x_row)
            y.append(int(r[1]))
            r = data_file.readline()
        x, y = np.array(x), np.array(y)
        x_train, x_test, y_train, y_test = cross_validation.train_test_split(x, y, test_size=cross_ratio, random_state=42)
        k_fold = 5
        kf = StratifiedKFold(y_train, n_folds=k_fold)
        results = []
        for c in [1]:
            for degree in [1]:
                cc = 0
                for train, test in kf:
                    cv_x_train, cv_x_test, cv_y_train, cv_y_test = x_train[train], x_train[test], y_train[train], y_train[test]
                    cf = svm.SVC(kernel='rbf', C=c, degree=degree, probability=True, random_state=42).fit(cv_x_train, cv_y_train)
                    predictions = cf.predict(cv_x_test)
                    pb = cf.predict_proba(cv_x_test)
                    comps = []
                    for a,b in zip(cv_y_test, predictions):
                        comps.append("%s-%s" % (a, b))
                    score = cf.score(cv_x_test, cv_y_test)
                    fpr, tpr, thresholds = roc_curve(cv_y_test, pb[:, 1])
                    roc_auc = auc(fpr, tpr)
                    cc += roc_auc
                cc = cc / float(k_fold)
                results.append((cc, c, degree))
        results.sort(reverse=True)
        best_c, best_d = results[0][1], results[0][2]
        cf = svm.SVC(kernel='rbf', probability=True, random_state=42).fit(x_train, y_train)
        y_score = cf.decision_function(x_test)
        fpr, tpr, thresholds = roc_curve(y_test, y_score, pos_label=1)
        auc_val = auc(fpr, tpr)
        f = open("%s/%s.results.tab" % (model_folder, fn), "at")
        row = ["svm", str(auc_val)]
        f.write("\t".join(row) + "\n")
        f.close()
        ncontrol, ncase = len(x[y==0]), len(x[y==1])
        save_roc(comps_id, "%s.svm" % fn, tpr, fpr, auc_val, ncontrol, ncase)

def predict_randomf(comps_id):
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    for fn, _, _ in types:
        print "%s:model:randomf.%s" % (comps_id, fn)
        data_file = open(os.path.join(model_folder, "%s.tab" % fn), "rt")
        header = data_file.readline().replace("\r", "").replace("\n", "").split("\t")
        r = data_file.readline()
        y = []
        x = []
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            x_row = []
            for el in r[2:]:
                x_row.append(float(el))
            x.append(x_row)
            y.append(int(r[1]))
            r = data_file.readline()
        x, y = np.array(x), np.array(y)
        x_train, x_test, y_train, y_test = cross_validation.train_test_split(x, y, test_size=cross_ratio, random_state=42)
        cf = ExtraTreesClassifier(n_estimators=250, random_state=42)
        cf.fit(x_train, y_train)

        importances = cf.feature_importances_
        std = np.std([tree.feature_importances_ for tree in cf.estimators_], axis=0)
        indices = np.argsort(importances)[::-1]
        print("Feature ranking:")
        for f in range(10):
            print("%d. feature %d %s (%f)" % (f + 1, indices[f], header[indices[f]+2], importances[indices[f]]))

        pb = cf.predict_proba(x_test)
        fpr, tpr, thresholds = roc_curve(y_test, pb[:, 1], pos_label=1)
        auc_val = auc(fpr, tpr)
        f = open("%s/%s.results.tab" % (model_folder, fn), "at")
        row = ["randomf", str(auc_val)]
        f.write("\t".join(row) + "\n")
        f.close()
        ncontrol, ncase = len(x[y==0]), len(x[y==1])
        save_roc(comps_id, "%s.randomf" % fn, tpr, fpr, auc_val, ncontrol, ncase)

        # plot feature ranking
        matplotlib.rcParams['axes.labelsize'] = 14
        matplotlib.rcParams['axes.titlesize'] = 14
        matplotlib.rcParams['xtick.labelsize'] = 12
        matplotlib.rcParams['ytick.labelsize'] = 12
        matplotlib.rcParams['legend.fontsize'] = 19
        matplotlib.rc('axes',edgecolor='gray')
        matplotlib.rcParams['axes.linewidth'] = 0.3
        matplotlib.rcParams['legend.frameon'] = 'False'
        pl.clf()
        a = pl.axes([0.1, 0.4, 0.85, 0.5])
        a.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        a.spines["top"].set_alpha(0.1)
        a.spines["bottom"].set_alpha(0.1)
        a.spines["left"].set_alpha(0.1)
        a.spines["right"].set_alpha(0.1)

        # version 1: plots variation inside forests
        #pl.bar(range(10), importances[indices][:10], color="g", yerr=std[indices][:10], align="center", edgecolor = "none")
        # version 2: no variation bar inside forests
        pl.bar(range(10), importances[indices][:10], color="lightgray", align="center", edgecolor = "none", width=0.5)
        pl.xticks(range(10), [header[indices[f]+2] for f in range(0, 10)], rotation='vertical')
        pl.xlim([-1, 10])
        pl.xlabel('Feature')
        pl.ylabel('Importance')
        pl.title('Feature importances')
        model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
        pl.savefig(os.path.join(model_folder, "%s_randomf_features.png" % fn))
        pl.savefig(os.path.join(model_folder, "%s_randomf_features.pdf" % fn))

def predict_tree(comps_id):
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    for fn in ["siteup_e", "siteup_r", "sitedown_e", "sitedown_r"]:
        print "%s:model:tree.%s" % (comps_id, fn)
        data_file = open(os.path.join(model_folder, "%s.tab" % fn), "rt")
        header = data_file.readline().replace("\r", "").replace("\n", "").split("\t")
        r = data_file.readline()
        y = []
        x = []
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            x_row = []
            for el in r[2:]:
                x_row.append(float(el))
            x.append(x_row)
            y.append(int(r[1]))
            r = data_file.readline()
        x, y = np.array(x), np.array(y)
        x_train, x_test, y_train, y_test = cross_validation.train_test_split(x, y, test_size=cross_ratio, random_state=42)
        cf = tree.DecisionTreeClassifier().fit(x_train, y_train)
        cf.predict
        pb = cf.predict_proba(x_test)
        fpr, tpr, thresholds = roc_curve(y_test, pb[:, 1])
        auc_val = auc(fpr, tpr)
        f = open("%s/%s.results.tab" % (model_folder, fn), "at")
        row = ["tree", str(auc_val)]
        f.write("\t".join(row) + "\n")
        f.close()
        #ncontrol, ncase = len(x_test[y_test==0]), len(x_test[y_test==1])
        ncontrol, ncase = len(x[y==0]), len(x[y==1])
        save_roc(comps_id, "%s.tree" % fn, tpr, fpr, auc_val, ncontrol, ncase)

def predict_linreg(comps_id, version="randomf"):
    print "%s:model:linreg" % comps_id
    from sklearn import datasets, linear_model
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    data_file = open(os.path.join(model_folder, "siteup.r.tab"), "rt")
    header = data_file.readline().replace("\r", "").replace("\n", "").split("\t")
    r = data_file.readline()
    y = []
    x = []
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        x_row = []
        for el in r[2:]:
            x_row.append(float(el))
        x.append(x_row)
        y.append(int(r[1]))
        r = data_file.readline()
    x, y = np.array(x), np.array(y)
    x_train, x_test, y_train, y_test = cross_validation.train_test_split(x, y, test_size=cross_ratio, random_state=42)
    cf = linear_model.LinearRegression().fit(x_train, y_train)
    print("Residual sum of squares: %.2f" % np.mean((cf.predict(x_test) - y_test) ** 2))
    print('Variance score: %.2f' % cf.score(x_test, y_test, ))

def plot_pca(comps_id):
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    for fn, _, _ in types:
        print "%s:model:pca.%s" % (comps_id, fn)
        data_file = open(os.path.join(model_folder, "%s.tab" % fn), "rt")
        header = data_file.readline().replace("\r", "").replace("\n", "").split("\t")
        r = data_file.readline()
        y = []
        x = []
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            x_row = []
            for el in r[2:]:
                x_row.append(float(el))
            x.append(x_row)
            y.append(int(r[1]))
            r = data_file.readline()
        x, y = np.array(x), np.array(y)
        pca = PCA(n_components=2)
        x_r = pca.fit(x).transform(x)
        # print len(x_r[y==1]) # this actually works
        pl.clf()
        a = pl.axes([0.1, 0.1, 0.85, 0.8])
        a.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        a.spines["top"].set_alpha(0.1)
        a.spines["bottom"].set_alpha(0.1)
        a.spines["left"].set_alpha(0.1)
        a.spines["right"].set_alpha(0.1)
        if fn.endswith(".e"):
            class_name = "enhanced"
        else:
            class_name = "repressed"
        for c, i, target_name in zip("rb", [0, 1], ["control", class_name]):
            pl.scatter(x_r[y == i, 0], x_r[y == i, 1], c=c, label=target_name, alpha=0.5, edgecolors='none')
        pl.legend()
        control, case = len(x[y==0]), len(x[y==1])
        pl.title('PCA %s : %s (#%s vs #%s)' % (fn, comps_id, case, control))
        pl.legend(loc="lower right")
        pl.savefig(os.path.join(model_folder, "pca.%s.png" % fn))
        pl.savefig(os.path.join(model_folder, "pca.%s.pdf" % fn))

def plot_lda(comps_id):
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    for fn in ["siteup_e", "siteup_r", "sitedown_e", "sitedown_r"]:
        print "%s:model:lda.%s" % (comps_id, fn)
        data_file = open(os.path.join(model_folder, "%s.tab" % fn), "rt")
        header = data_file.readline().replace("\r", "").replace("\n", "").split("\t")
        r = data_file.readline()
        y = []
        x = []
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            x_row = []
            for el in r[2:]:
                x_row.append(float(el))
            x.append(x_row)
            y.append(int(r[1]))
            r = data_file.readline()
        x, y = np.array(x), np.array(y)
        lda = LDA(n_components=2)
        x_r = lda.fit(x, y).transform(x)
        # print len(x_r[y==1]) # this actually works
        pl.clf()
        a = pl.axes([0.1, 0.1, 0.85, 0.8])
        a.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        a.spines["top"].set_alpha(0.1)
        a.spines["bottom"].set_alpha(0.1)
        a.spines["left"].set_alpha(0.1)
        a.spines["right"].set_alpha(0.1)
        if fn.endswith(".e"):
            class_name = "enhanced"
        else:
            class_name = "repressed"
        for c, i, target_name in zip("rb", [0, 1], ["control", class_name]):
            pl.scatter(x_r[y == i, 0], x_r[y == i, 1], c=c, label=target_name, alpha=0.5, edgecolors='none')
        pl.legend()
        control, case = len(x[y==0]), len(x[y==1])
        pl.title('LDA %s : %s (#%s vs #%s)' % (fn, comps_id, case, control))
        pl.legend(loc="lower right")
        pl.savefig(os.path.join(model_folder, "lda.%s.png" % fn))
        pl.savefig(os.path.join(model_folder, "lda.%s.pdf" % fn))

def plot_roc(comps_id):
    print "%s:model:plot_roc" % comps_id
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    pickle_file = os.path.join(model_folder, "roc.pickle")
    auc = pickle.load(open(pickle_file, "rb"))

    # styling
    matplotlib.rcParams['axes.labelsize'] = 17
    matplotlib.rcParams['axes.titlesize'] = 17
    matplotlib.rcParams['xtick.labelsize'] = 14
    matplotlib.rcParams['ytick.labelsize'] = 14
    matplotlib.rcParams['legend.fontsize'] = 14
    matplotlib.rc('axes',edgecolor='gray')
    matplotlib.rcParams['axes.linewidth'] = 0.3
    matplotlib.rcParams['legend.frameon'] = 'False'

    for fn, _, _ in types:
        pl.clf()
        a = pl.axes([0.1, 0.1, 0.85, 0.8])
        a.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        a.spines["top"].set_alpha(0.1)
        a.spines["bottom"].set_alpha(0.1)
        a.spines["left"].set_alpha(0.1)
        a.spines["right"].set_alpha(0.1)
        control, case = 0, 0
        for version, (x, y, auc_val, ncontrol, ncase) in auc.items():
            if version.startswith(fn):
                if control==0:
                    control = ncontrol
                else:
                    assert(control==ncontrol)
                if case==0:
                    case = ncase
                else:
                    assert(case==ncase)
                pl.plot(x, y, label='%s : AUC = %0.3f' % (version.split(".")[-1], auc_val))
        pl.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        pl.xlim([0.0, 1.0])
        pl.ylim([0.0, 1.0])
        pl.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0", "0.2", "0.4", "0.6", "0.8", "1.0"])
        pl.yticks([0.2, 0.4, 0.6, 0.8, 1], ["0.2", "0.4", "0.6", "0.8", "1.0"])

        pl.xlabel('False Positive Rate')
        pl.ylabel('True Positive Rate')
        pl.title('ROC: %s : %s (#%s vs #%s)' % (fn, comps_id, case, control))
        pl.legend(loc="lower right")
        pl.savefig(os.path.join(model_folder, "roc.%s.png" % fn))
        pl.savefig(os.path.join(model_folder, "roc.%s.pdf" % fn))

def save_roc(comps_id, version, tpr, fpr, auc_val, ncontrol, ncase):
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    pickle_file = os.path.join(model_folder, "roc.pickle")
    if os.path.exists(pickle_file):
        auc = pickle.load(open(pickle_file, "rb"))
    else:
        auc = {}
    auc[version] = (fpr, tpr, auc_val, ncontrol, ncase)
    pickle.dump(auc, open(pickle_file, "wb"))

# http://sebastianraschka.com/Articles/2014_python_lda.html
def plot_hist(comps_id):
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    for fn in ["siteup.r"]: #, "siteup.e", "sitedown.r", "sitedown.e"]:
        data_file = open(os.path.join(model_folder, "%s.tab" % fn), "rt")
        header = data_file.readline().replace("\r", "").replace("\n", "").split("\t")
        r = data_file.readline()
        y = []
        x = []
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            x_row = []
            for el in r[2:]:
                x_row.append(float(el))
            x.append(x_row)
            y.append(int(r[1]))
            r = data_file.readline()
        x, y = np.array(x), np.array(y)

    at = header.index("mTGTG")-2
    min_b = math.floor(np.min(x[:,at]))
    max_b = math.ceil(np.max(x[:,at]))
    bins = np.linspace(min_b, max_b, 25)

    pl.clf()
    a = pl.axes([0.1, 0.1, 0.85, 0.8])
    a.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    a.spines["top"].set_alpha(0.1)
    a.spines["bottom"].set_alpha(0.1)
    a.spines["left"].set_alpha(0.1)
    a.spines["right"].set_alpha(0.1)
    # plottling the histograms
    for lab, lab_name, col in zip([0, 1], ("control", "regulated"), ('blue', 'red')):
        a.hist(x[y==lab, at],
                   color=col,
                   label='class %s' % lab_name,
                   bins=bins,
                   alpha=0.5)
    ylims = a.get_ylim()

    # plot annotation
    leg = a.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    a.set_ylim([0, max(ylims)+2])
    a.set_xlabel("feature1")
    a.set_title("hist")

    # hide axis ticks
    a.tick_params(axis="both", which="both", bottom="off", top="off",
            labelbottom="on", left="off", right="off", labelleft="on")
    pl.savefig(os.path.join(model_folder, "hist.%s.png" % fn))
    pl.savefig(os.path.join(model_folder, "hist.%s.pdf" % fn))

def write_index(comps_id):
    print "%s:model:write_index" % comps_id
    model_folder = os.path.join(apa.path.comps_folder, comps_id, "model")
    pickle_file = os.path.join(model_folder, "roc.pickle")
    auc = pickle.load(open(pickle_file, "rb"))
    fn = os.path.join(model_folder, "model.html")
    f = open(fn, "wt")
    f.write("<html>\n")
    head = """<head>
    <script type="text/javascript" src="../../software/js/jquery-1.8.0.min.js"></script>
    <script type="text/javascript" src="../../software/js/jquery-ui-1.8.23.custom.min.js"></script>
    <link rel="stylesheet" type="text/css" href="../../software/tooltipster-master/css/tooltipster.css" />
    <script type="text/javascript" src="../../software/tooltipster-master/js/jquery.tooltipster.min.js"></script>
    </head>"""
    f.write(head+"\n")
    f.write("<body>\n")
    f.write("<center>")

    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 16px;'><tr>")
    for fn, _, _ in types:
        f.write("<td align=center>PCA %s</td>" % fn)
    f.write("</tr>\n")
    #for fn in ["siteup.r", "siteup.e", "sitedown.r", "sitedown.e"]:
    for fn, _, _ in types:
        f.write("<td align=right valign=center width=30px>")
        f.write("<a href=%s><img src=%s width=400px></a>" % ("pca.%s.png" % fn, "pca.%s.png" % fn))
        f.write("</td>")
    f.write("</tr>")
    f.write("\n")
    f.write("</table>")

    f.write("<br>")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 16px;'><tr>")
    for fn, _, _ in types:
        text = ""
        for t in models:
            (fpr, tpr, auc_val, ncontrol, ncase) = auc["%s.%s" % (fn, t)]
            text += "<br>AUC %s = %.3f" % (t, auc_val)
        f.write("<td align=center>ROC %s<br>%s</td>" % (fn, text))
    f.write("</tr>\n")
    for fn, _, _ in types:
        f.write("<td align=right valign=center width=30px>")
        f.write("<a href=%s><img src=%s width=400px></a>" % ("roc.%s.png" % fn, "roc.%s.png" % fn))
        f.write("</td>")
    f.write("</tr>")
    f.write("\n")
    f.write("</table>")

    f.write("<br>")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 16px;'><tr>")
    for fn, _, _ in types:
        f.write("<td align=center>RandomF %s</td>" % fn)
    f.write("</tr>\n")
    for fn, _, _ in types:
        f.write("<td align=right valign=center width=30px>")
        f.write("<a href=%s><img src=%s width=400px></a>" % ("%s_randomf_features.png" % fn, "%s_randomf_features.png" % fn))
        f.write("</td>")
    f.write("</tr>")
    f.write("\n")
    f.write("</table>")

    f.write("<br>")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 16px;'><tr><td align=center>siteup e vs c</td><td align=center width=15px>siteup r vs c</td></tr>\n")
    f.write("<tr>")
    for fn in ["siteup_e", "siteup_r"]:
        f.write("<td align=right valign=center width=30px>")
        f.write("<a href=%s><img src=%s width=400px></a>" % ("%s_ntdist.png" % fn, "%s_ntdist.png" % fn))
        f.write("</td>")
    f.write("</tr>")
    f.write("\n")
    f.write("</table>")

    f.write("<br>")
    f.write("<table style='border-collapse: collapse; border-spacing: 0px; font-size: 16px;'><tr><td align=center>sitedown e vs c</td><td align=center width=15px>sitedown r vs c</td></tr>\n")
    f.write("<tr>")
    for fn in ["sitedown_e", "sitedown_r"]:
        f.write("<td align=right valign=center width=30px>")
        f.write("<a href=%s><img src=%s width=400px></a>" % ("%s_ntdist.png" % fn, "%s_ntdist.png" % fn))
        f.write("</td>")
    f.write("</tr>")
    f.write("\n")
    f.write("</table>")

    f.write("</body>")
    f.write("</html>\n")
    f.close()
