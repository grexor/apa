import pybio
d = pybio.data.Bedgraph("data.bed") # load initial data
d.load("data2.bed")                 # load additional data
d.get_region("chr1", "+", 100)      # get value at loci +chr1:100
