# main apa module

import sys
import os
import pybio

import apa
import config
import path
import annotation
import extract
import map
import bed
import polya
import comps
import rnamap
import model
import analysis

apa.path.init() # inicialize paths
apa.annotation.init() # read annotations
