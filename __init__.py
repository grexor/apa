"""
apa
====

Research platform for bioinformatics analysis of NGS alternative polyadenylation data.
"""
import sys
import os
import pybio

import apa.config
import apa.analysis
import apa.path
import apa.annotation
import apa.extract
import apa.map
import apa.bed
import apa.polya
import apa.comps
import apa.rnamap
#import model # to-do
import apa.motifs

version = "1.3"

apa.path.init() # inicialize paths
apa.config.init()
apa.annotation.init() # read annotations

def get_read_len(lib_id, exp_id):
    f = pybio.data.Fastq(apa.path.map_fastq_file(lib_id, exp_id))
    f.read()
    return len(f.sequence)
