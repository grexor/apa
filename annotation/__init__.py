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
        apa.annotation.libs[lib_id] = apa.annotation.read(filename)

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
        self.dcode_len = 0 # demultiplex barcode length (first_5), the rest of the read is random barcode

def read(lib_id):
    """
    | Read library annotation, including experiments.

    | This is read from the following tab delimited file:
    |   :green:`data_folder/lib_id/annotation.tab`

    The structure of the annotation.tab file is the same for all libraries:

    ====== ======= === ====== ======= ======= ======
    exp_id method  rep tissue cond    species map_to
    ====== ======= === ====== ======= ======= ======
    1      lex_fwd 1   HeLa   control Hs      hg19
    2      lex_fwd 2   HeLa   control Hs      hg19
    ...
    ====== ======= === ====== ======= ======= ======

    The description of TAB columns:

    ============= ===========
    Column        Description
    ============= ===========
    exp_id        starts with 1 for each library
    method        3' end sequencing protocol
    rep           replicate number
    tissue        description of tissue the sample was taken from
    cond          description of the experimental conditions
    species       sample species
    map_to        genome assembly for mapping experiment's sequences
    ============= ===========
    """
    lib = Library(lib_id)
    data = {}
    filename = os.path.join(apa.path.data_folder, lib_id, "annotation.tab")
    f = open(filename, "rt")
    lib.fastq_files = f.readline().replace("\r", "").replace("\n", "").replace("\t", "").split("&")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()

    dcode_len = None
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        try:
            dcode = data["first_5"].split(",")[0].split("_")[0] # take exact demulti code
            # all experiments in the library have the same first_5 code length and random barcode length
            if dcode_len==None:
                dcode_len = len(dcode)
            else:
                assert(dcode_len == len(dcode))
        except:
            dcode = None
            dcode_len = None
        try:
            rcode = data["first_5"].split(",")[1] # take exact random code
        except:
            rcode = None

        data["dcode"] = dcode
        lib.experiments[int(data["exp_id"])] = data
        r = f.readline()
    lib.dcode_len = dcode_len
    return lib
