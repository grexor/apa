import apa
import os
import pybio
import glob
import struct
import fcntl
import random
import errno
import time

def init():
    apa.annotation.libs = {}
    files = glob.glob(os.path.join(apa.path.data_folder, "*"))
    for filename in files:
        if not os.path.isdir(filename):
            continue
        lib_id = filename.split("/")[-1]
        annotation_tab = os.path.join(apa.path.data_folder, lib_id, "annotation.tab")
        if os.path.exists(annotation_tab):
            apa.annotation.libs[lib_id] = apa.annotation.Library(lib_id)

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
        self.owner = []
        self.access = []
        self.name = ""
        self.notes = ""
        self.tags = ""
        self.genome = ""
        self.map_to = ""
        self.method = ""
        self.public_only = [];
        self.seq_type = "single" # sequencing type: single, paired
        self.columns = [("Tissue", "tissue"), ("Condition", "condition"), ("Replicate", "replicate"), ("Upload Filename_R1", "upload_filename_r1"), ("Upload Filename_R2", "upload_filename_r2")]
        self.columns_display = [("Tissue", "tissue"), ("Condition", "condition"), ("Replicate", "replicate"), ("Upload Filename_R1", "upload_filename_r1"), ("Upload Filename_R2", "upload_filename_r2")]
        self.authors = []
        self.status = ""
        if lib_id!=None:
            self.read_lib(lib_id)

    def read_lib(self, lib_id):
        data = {}
        filename = os.path.join(apa.path.data_folder, lib_id, "annotation.tab")
        if not os.path.exists(filename):
            return
        while True:
            f = open(filename, "rt")
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except IOError as e:
                # raise on unrelated IOErrors
                if e.errno != errno.EAGAIN:
                    raise
                else:
                    time.sleep(random.random())
        self.fastq_files = f.readline().replace("\r", "").replace("\n", "").replace("\t", "").split("&")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            self.experiments[int(data["exp_id"])] = data
            r = f.readline()
        fcntl.flock(f, fcntl.LOCK_UN)
        f.close()

        filename = os.path.join(apa.path.data_folder, lib_id, "%s.config" % lib_id)
        if os.path.exists(filename):
            while True:
                f = open(filename, "rt")
                try:
                    fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except IOError as e:
                    # raise on unrelated IOErrors
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(random.random())
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
                if r[0].startswith("status:"):
                    self.status = str(r[0].split("status:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("access:"):
                    self.access = str(r[0].split("access:")[1]).split(",")
                    r = f.readline()
                    continue
                if r[0].startswith("genome:"):
                    self.genome = str(r[0].split("genome:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("map_to:"):
                    self.map_to = str(r[0].split("map_to:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("method:"):
                    self.method = str(r[0].split("method:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("seq_type:"):
                    self.seq_type = str(r[0].split("seq_type:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("authors:"):
                    self.authors = str(r[0].split("authors:")[1]).split(",")
                    r = f.readline()
                    continue
                if r[0].startswith("public_only:"):
                    self.public_only = str(r[0].split("public_only:")[1]).split(",")
                    r = f.readline()
                    continue
                if r[0].startswith("owner:"):
                    self.owner = str(r[0].split("owner:")[1]).split(",")
                    r = f.readline()
                    continue
                if r[0].startswith("name:"):
                    self.name = str(r[0].split("name:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("notes:"):
                    self.notes = str(r[0].split("notes:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("tags:"):
                    self.tags = str(r[0].split("tags:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("columns:"):
                    self.columns = eval(r[0].split("columns:")[1])
                    r = f.readline()
                    continue
                if r[0].startswith("columns_display:"):
                    self.columns_display = eval(r[0].split("columns_display:")[1])
                    r = f.readline()
                    continue
                r = f.readline()
            fcntl.flock(f, fcntl.LOCK_UN)
            f.close()
            if self.map_to=="":
                self.map_to = self.genome

    def add_status(self, status_name=None):
        if self.status!="":
            status = set(self.status.split(","))
        else:
            status = set()
        if status!=None:
            status.add(status_name)
        self.status = ",".join(list(status))

    def remove_status(self, status_name):
        status = self.status.split(",")
        status = filter(lambda el: el != status_name, status)
        self.status = ",".join(status)

    def save(self):
        filename = os.path.join(apa.path.data_folder, self.lib_id, "%s.config" % self.lib_id)
        while True:
            f = open(filename, "wt")
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except IOError as e:
                # raise on unrelated IOErrors
                if e.errno != errno.EAGAIN:
                    raise
                else:
                    time.sleep(random.random()*2+1)
        f.write("access:" + ",".join(self.access) + "\n")
        f.write("name:%s" % (self.name) + "\n")
        f.write("notes:%s" % (self.notes) + "\n")
        f.write("public_only:" + ",".join(self.public_only) + "\n")
        f.write("owner:" + ",".join(self.owner) + "\n")
        f.write("method:" + self.method + "\n")
        f.write("genome:" + self.genome + "\n")
        f.write("map_to:" + self.map_to + "\n")
        f.write("seq_type:" + self.seq_type + "\n")
        f.write("tags:" + self.tags + "\n")
        f.write("status:" + self.status + "\n")
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
        try:
            fcntl.flock(f, fcntl.LOCK_UN)
            f.close()
        except:
            pass

        filename = os.path.join(apa.path.data_folder, self.lib_id, "annotation.tab")
        while True:
            f = open(filename, "wt")
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except IOError as e:
                # raise on unrelated IOErrors
                if e.errno != errno.EAGAIN:
                    raise
                else:
                    time.sleep(random.random()*2+1)
        f.write("%s\n" % self.lib_id)
        columns = ["exp_id", "species", "map_to", "method"]
        for (cname, cid) in self.columns:
            columns.append(cid)
        f.write("%s\n" % ("\t".join(columns)))
        exp_ids = list(self.experiments.keys())
        exp_ids = sorted(exp_ids, key=lambda x: (int(x))) # sort, but keep "strings" in case they are string
        for exp_id in exp_ids:
            row = [exp_id, self.genome, self.map_to, self.method]
            for cid in columns[4:]:
                try:
                    row.append(self.experiments[exp_id][cid])
                except:
                    row.append("")
            f.write("%s\n" % ("\t".join([str(el) for el in row])))
        try:
            fcntl.flock(f, fcntl.LOCK_UN)
            f.close()
        except:
            pass
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
