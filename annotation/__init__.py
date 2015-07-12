import apa
import os
import pybio
import glob
import struct

def init():
    apa.annotation.libs = {}
    files = glob.glob(os.path.join(apa.path.data_folder, "*"))
    for filename in files:
        library_id = filename.split("/")[-1]
        apa.annotation.libs[library_id] = apa.annotation.read(filename)

def aremoved(library_id, read_id):
    filename = os.path.join(apa.path.data_folder, library_id, "%s.aremoved.bin" % (library_id))
    if not os.path.exists(filename):
        return None
    aremoved_file = open(filename, "rb")
    aremoved_file.seek(read_id)
    data = aremoved_file.read(1)
    if len(data)==0:
        return 0
    return struct.unpack("B", data)[0]

def rndcode(library_id, read_id):
    rndcode_file = open(os.path.join(apa.path.data_folder, library_id, "%s.rnd.bin" % (library_id)), "rb")
    rndcode_file.seek(read_id*4)
    data = rndcode_file.read(4)
    if len(data)==0:
        return 0
    return struct.unpack("I", data)[0]

class Library:  
    def __init__(self, library_id):
        self.library_id = library_id
        self.experiments = {}
        self.fastq_files = {}
        self.dcode_len = 0 # demultiplex barcode length (first_5), the rest of the read is random barcode

def read(library_id):
    lib = Library(library_id)
    data = {}
    filename = os.path.join(apa.path.data_folder, library_id, "annotation.tab")
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