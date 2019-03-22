import apa
import os
import pybio
import glob
import struct

def init():
    apa.annotation.libs = {}
    files = glob.glob(os.path.join(apa.path.data_folder, "*"))
    for filename in files:
        if not os.path.isdir(filename):
            continue
        lib_id = filename.split("/")[-1]
        annotation_tab = os.path.join(apa.path.data_folder, lib_id, "annotation.tab")
        if os.path.exists(annotation_tab):
            apa.annotation.libs[lib_id] = apa.annotation.read(lib_id)

def aremoved(lib_id, read_id):
    filename = os.path.join(apa.path.data_folder, lib_id, "%s.aremoved.bin" % (lib_if))
    if not os.path.exists(filename):
        return None
    aremoved_file = open(filename, "rb")
    aremoved_file.seek(read_id)
    data = aremoved_file.read(1)
    if len(data)==0:
        return 0
    return struct.unpack("B", data)[0]

def rndcode(lib_id, read_id):
    rndcode_file = open(os.path.join(apa.path.data_folder, lib_id, "%s.rnd.bin" % (lib_id)), "rb")
    rndcode_file.seek(read_id*4)
    data = rndcode_file.read(4)
    if len(data)==0:
        return 0
    return struct.unpack("I", data)[0]

class Library:
    def __init__(self, lib_id):
        self.lib_id = lib_id
        self.experiments = {}
        self.fastq_files = {}
        self.owner = [];
        self.access = [];
        self.name = "";
        self.notes = "";
        self.genome = "";
        self.method = ""
        self.public_only = [];
        self.seq_type = "single" # sequencing type: single, paired
        self.columns = [("Tissue", "tissue"), ("Condition", "condition"), ("Replicate", "replicate"), ("Upload Filename_R1", "upload_filename_r1"), ("Upload Filename_R2", "upload_filename_r2")]
        self.columns_display = [("Tissue", "tissue"), ("Condition", "condition"), ("Replicate", "replicate"), ("Upload Filename_R1", "upload_filename_r1"), ("Upload Filename_R2", "upload_filename_r2")]
        self.authors = []

    def save(self):
        filename = os.path.join(apa.path.data_folder, self.lib_id, "%s.config" % self.lib_id)
        f = open(filename, "wt")
        f.write("access:" + ",".join(self.access) + "\n")
        f.write("name:%s" % (self.name) + "\n")
        f.write("notes:%s" % (self.notes) + "\n")
        f.write("public_only:" + ",".join(self.public_only) + "\n")
        f.write("owner:" + ",".join(self.owner) + "\n")
        f.write("method:" + self.method + "\n")
        f.write("genome:" + self.genome + "\n")
        f.write("seq_type:" + self.seq_type + "\n")
        columns = []
        for c in self.columns:
            if self.seq_type=="single" and c[1]=="upload_filename_r2":
                continue
            columns.append(c)
        columns_display = []
        for c in self.columns_display:
            if self.seq_type=="single" and c[1]=="upload_filename_r2":
                continue
            columns_display.append(c)
        f.write("columns:%s" % (str(columns)) + "\n")
        f.write("columns_display:%s" % (str(columns_display)) + "\n")
        filename = os.path.join(apa.path.data_folder, self.lib_id, "annotation.tab")
        f = open(filename, "wt")
        f.write("%s\n" % self.lib_id)
        columns = ["exp_id", "species", "map_to", "method"]
        for (cname, cid) in self.columns:
            columns.append(cid)
        f.write("%s\n" % ("\t".join(columns)))
        exp_ids = self.experiments.keys()
        exp_ids = sorted(exp_ids, key=lambda x: (int(x))) # sort, but keep "strings" in case they are string
        for exp_id in exp_ids:
            row = [exp_id, self.genome, self.genome, self.method]
            for cid in columns[4:]:
                try:
                    row.append(self.experiments[exp_id][cid])
                except:
                    row.append("")
            f.write("%s\n" % ("\t".join([str(el) for el in row])))
        f.close()
        return True

    def add_experiment(self, data):
        exp_id = len(self.experiments)+1
        self.experiments[exp_id] = data
        return exp_id

    def add_empty_experiment(self, filename_R1="", filename_R2=""):
        exp_id = len(self.experiments)+1
        data = {"species":"", "map_to":""}
        for _, cid in self.columns:
            data[cid] = ""
        data["upload_filename_r1"] = filename_R1
        data["upload_filename_r2"] = filename_R2
        self.experiments[exp_id] = data
        return exp_id

    def edit_experiment(self, exp_id, data):
        self.experiments[exp_id] = data

def count_ownership(email):
    num_libs = 0
    num_experiments = 0
    for exp_id, data in apa.annotation.libs.items():
        if email in data.owner:
            num_libs += 1
            num_experiments += len(data.experiments)
    return num_libs, num_experiments

def read(lib_id):
    lib = Library(lib_id)
    data = {}
    filename = os.path.join(apa.path.data_folder, lib_id, "annotation.tab")
    f = open(filename, "rt")
    lib.fastq_files = f.readline().replace("\r", "").replace("\n", "").replace("\t", "").split("&")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        lib.experiments[int(data["exp_id"])] = data
        r = f.readline()

    filename = os.path.join(apa.path.data_folder, lib_id, "%s.config" % lib_id)
    if os.path.exists(filename):
        f = open(filename, "rt")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "")
            r = r.split(" #")[0] # remove comments
            r = r.rstrip() # remove whitespace characters from end of string
            r = r.split("\t")
            if r==[""]:
                r = f.readline()
                continue
            if r[0].startswith("#"):
                r = f.readline()
                continue
            if r[0].startswith("access:"):
                lib.access = str(r[0].split("access:")[1]).split(",")
                r = f.readline()
                continue
            if r[0].startswith("genome:"):
                lib.genome = str(r[0].split("genome:")[1])
                r = f.readline()
                continue
            if r[0].startswith("method:"):
                lib.method = str(r[0].split("method:")[1])
                r = f.readline()
                continue
            if r[0].startswith("seq_type:"):
                lib.seq_type = str(r[0].split("seq_type:")[1])
                r = f.readline()
                continue
            if r[0].startswith("authors:"):
                lib.authors = str(r[0].split("authors:")[1]).split(",")
                r = f.readline()
                continue
            if r[0].startswith("public_only:"):
                lib.public_only = str(r[0].split("public_only:")[1]).split(",")
                r = f.readline()
                continue
            if r[0].startswith("owner:"):
                lib.owner = str(r[0].split("owner:")[1]).split(",")
                r = f.readline()
                continue
            if r[0].startswith("name:"):
                lib.name = str(r[0].split("name:")[1])
                r = f.readline()
                continue
            if r[0].startswith("notes:"):
                lib.notes = str(r[0].split("notes:")[1])
                r = f.readline()
                continue
            if r[0].startswith("columns:"):
                lib.columns = eval(r[0].split("columns:")[1])
                r = f.readline()
                continue
            if r[0].startswith("columns_display:"):
                lib.columns_display = eval(r[0].split("columns_display:")[1])
                r = f.readline()
                continue
            r = f.readline()
        f.close()
    return lib
