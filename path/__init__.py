import apa
import os
import glob

data_folder = "data.apa"
polya_folder = "data.polya"
comps_folder = "data.comps"
iCLIP_folder = "data.iCLIP"
rnamotifs_folder = "/home/gregor/rnamotifs2"

def init():
    apa.path.root_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", ".."))
    apa.path.data_folder = os.path.join(root_folder, apa.path.data_folder)
    apa.path.comps_folder = os.path.join(root_folder, apa.path.comps_folder)
    apa.path.iCLIP_folder = os.path.join(root_folder, apa.path.iCLIP_folder)
    apa.path.polya_folder = os.path.join(root_folder, apa.path.polya_folder)

def r_filename(lib_id, exp_id, map_id=1):
    """
    Returns constructed path to :ref:`R bedGraph file <r_bedgraph_method>` from lib_id, exp_id and map_id:

    .. code-block:: bash

        ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/lib_id_e${exp_id}_m${map_id}.raw.bg
    """
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.raw.bed.gz" % (lib_id, exp_id, map_id))

def bam_filename(lib_id, exp_id, map_id=1):
    """
    Returns constructed path to :ref:`R bedGraph file <r_bedgraph_method>` from lib_id, exp_id and map_id:

    .. code-block:: bash

        ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/lib_id_e${exp_id}_m${map_id}.bam
    """
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.bam" % (lib_id, exp_id, map_id))

def t_filename(lib_id, exp_id, map_id=1):
    """
    Returns constructed path to :ref:`tail bedGraph file <t_bedgraph_method>` from lib_id, exp_id and map_id:

    .. code-block:: bash

        ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/lib_id_e${exp_id}_m${map_id}.tail.bed.gz
    """
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s.tail.bed.gz" % (lib_id, exp_id, map_id))

def e_filename(lib_id, exp_id, map_id=1, poly_id=""):
    """
    Returns constructed path to :ref:`expression bedGraph file <e_bedgraph_method>` from lib_id, exp_id, poly_id and map_id:

    .. code-block:: bash

        ${data_folder}/${lib_id}/e${exp_id}/m${map_id}/lib_id_e${exp_id}_m${map_id}_db-${poly_id}.exp.bed.gz
    """
    if poly_id=="":
        poly_id = lib_id
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s_db-%s.exp.bed.gz" % (lib_id, exp_id, map_id, poly_id))

def e_filename_norm(lib_id, exp_id, map_id=1, poly_id=""):
    if poly_id=="":
        poly_id = lib_id
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s_e%s_m%s_db-%s.exp_norm.bed.gz" % (lib_id, exp_id, map_id, poly_id))

def lock_filename(lib_id, exp_id, lock="bed", map_id=1):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s" % map_id, "%s.lock" % lock)

# poly_type = "strong", "weak", "less", "noclass"
def polyadb_filename(poly_id, poly_type=None, filetype="bed"):
    if filetype=="bed":
        if poly_type!=None:
            return os.path.join(apa.path.polya_folder, "%s_%s.bed.gz" % (poly_id, poly_type))
        if poly_type==None or poly_type=="all":
            return os.path.join(apa.path.polya_folder, "%s.bed.gz" % poly_id)
    if filetype=="pas":
        return os.path.join(apa.path.polya_folder, "%s_pas.bed.gz" % poly_id)
    if filetype=="temp":
        return os.path.join(apa.path.polya_folder, "%s.temp.gz" % poly_id)
    if filetype=="complete":
        return os.path.join(apa.path.polya_folder, "%s_complete.tab.gz" % poly_id)
    if filetype=="tab":
        if poly_type!=None:
            return os.path.join(apa.path.polya_folder, "%s_%s.tab.gz" % (poly_id, poly_type))
        else:
            return os.path.join(apa.path.polya_folder, "%s.tab.gz" % poly_id)
    if filetype=="polyar_pdf":
        return os.path.join(apa.path.polya_folder, "%s.pdf" % poly_id)
    return os.path.join(apa.path.polya_folder, "%s.%s" % (poly_id, filetype))

def polyadb_ann_filename(species):
    return os.path.join(apa.path.polya_folder, "polyadb.%s.tab.gz" % species)

def lib_folder(lib_id):
    return os.path.join(apa.path.data_folder, lib_id)

def map_folder(lib_id, exp_id, map_id=1, append=""):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "m%s%s" % (map_id, append))

def map_fastq_file(lib_id, exp_id, append=""):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "%s_e%s%s.fastq.bz2" % (lib_id, exp_id, append))

def map_fastq_file_raw(lib_id, exp_id):
    return os.path.join(apa.path.data_folder, lib_id, "e%s" % exp_id, "%s_e%s_raw.fastq.bz2" % (lib_id, exp_id))

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
