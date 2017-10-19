"""
apa
====

Research platform for bioinformatics analysis of NGS alternative polyadenylation data.
"""
import sys
import os
import pybio

import apa
import config
import analysis
import path
import annotation
import extract
import map
import bed
import polya
import comps
import rnamap
#import model # update this module in version 2
import motifs
import warnings

version = "1.1"

apa.path.init() # inicialize paths
apa.annotation.init() # read annotations

warnings.simplefilter('ignore')

def get_read_len(lib_id, exp_id):
    f = pybio.data.Fastq(apa.path.map_fastq_file(lib_id, exp_id))
    f.read()
    return len(f.sequence)
