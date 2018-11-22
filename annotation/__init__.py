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
        self.columns = [("Tissue", "tissue"), ("Condition", "condition"), ("Replicate", "replicate")]
        self.authors = []

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
            r = f.readline()
        f.close()
    return lib

def save(library):
    filename = os.path.join(apa.path.data_folder, library.lib_id, "%s.config" % library.lib_id)
    f = open(filename, "wt")
    f.write("access:" + ",".join(library.access) + "\n")
    f.write("name:%s" % (library.name) + "\n")
    f.write("notes:%s" % (library.notes) + "\n")
    f.write("public_only:" + ",".join(library.public_only) + "\n")
    f.write("owner:" + ",".join(library.owner) + "\n")
    f.write("columns:%s" % (str(library.columns)) + "\n")
    filename = os.path.join(apa.path.data_folder, library.lib_id, "annotation.tab")
    if not os.path.exists(filename):
        open(filename, "wt").close()
    return True
