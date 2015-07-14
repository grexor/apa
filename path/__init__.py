import apa
import os
import glob

def init():
    apa.path.root_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", ".."))
    apa.path.data_folder = os.path.join(root_folder, apa.config.data_folder)
    apa.path.comps_folder = os.path.join(root_folder, apa.config.comps_folder)
    apa.path.iCLIP_folder = os.path.join(root_folder, apa.config.iCLIP_folder)
    apa.path.polya_folder = os.path.join(root_folder, apa.config.polya_folder)

def t_filename(lib_id, exp_id, map_id=1):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.T.bed" % (lib_id, exp_id, map_id))

def r_filename(lib_id, exp_id, map_id=1):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.R.bed" % (lib_id, exp_id, map_id))

def lock_filename(lib_id, exp_id, lock="bed", map_id=1):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s.lock" % lock)

def e_filename(lib_id, exp_id, map_id=1, polyid=None, filetype=None):
    if polyid==None:
        polyid = apa.annotation.libs[lib_id].experiments[exp_id]["map_to"]
    if filetype!=None:
        return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s_%s.E.%s.bed" % (lib_id, exp_id, map_id, polyid, filetype))
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s_%s.E.bed" % (lib_id, exp_id, map_id, polyid))

def polyadb_filename(species, filetype="bed"):
    if filetype=="bed":
        return os.path.join(apa.path.polya_folder, "%s.bed" % species)
    if filetype=="temp":
        return os.path.join(apa.path.polya_folder, "%s.temp" % species)
    if filetype=="complete":
        return os.path.join(apa.path.polya_folder, "%s_complete.tab" % species)
    if filetype=="tab":
        return os.path.join(apa.path.polya_folder, "%s.tab" % species)
    if filetype=="class_hist":
        return os.path.join(apa.path.polya_folder, "%s_class_hist" % species)
    return os.path.join(apa.path.polya_folder, "%s.%s" % (species, filetype))

def polyadb_ann_filename(species):
    return os.path.join(apa.path.polya_folder, "polyadb.%s.tab" % species)

def map_folder(lib_id, exp_id, map_id=1):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id)

def map_fastq_file(lib_id, exp_id):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "%s_e%s.fastq.gz" % (lib_id, exp_id))

# comps
def comps_config_filename(comps_id):
    return os.path.join(apa.path.comps_folder, comps_id, "%s.config" % comps_id)

def comps_filename(comps_id, filetype):
    return os.path.join(apa.path.comps_folder, comps_id, "%s.%s" % (comps_id, filetype))

def comps_expression_filename(comps_id, filetype="genes"):
    return os.path.join(apa.path.comps_folder, comps_id, "%s.expression_%s.tab" % (comps_id, filetype))

def data_expression_genes(lib_id):
    return os.path.join(apa.path.data_folder, lib_id, "%s.genes.tab" % lib_id)

def data_expression_sites(lib_id):
    return os.path.join(apa.path.data_folder, lib_id, "%s.sites.tab" % lib_id)
