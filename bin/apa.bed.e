#!/usr/bin/python
import apa
import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-lib_id', action="store", dest="lib_id")
parser.add_argument('-poly_id', action="store", dest="poly_id")
parser.add_argument('-force', action="store_true", default=False)
parser.add_argument('-genome', action="store", default=None)
args = parser.parse_args()

apa.bed.process_e_lib(args.lib_id, polyid=args.poly_id, force=args.force, genome=args.genome)
