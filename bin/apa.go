#!/usr/bin/python3
import apa
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-comps_id', action="store", dest="comps_id", default=None)
args = parser.parse_args()

apa.analysis.go_enrichment(args.comps_id)
