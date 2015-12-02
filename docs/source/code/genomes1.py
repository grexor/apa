import pybio
pybio.genomes.load("hg19")                                # load hg19 genome annotation
pybio.genomes.annotate_position("hg19", 1, "+", 122434)   # (u'ENSG00000186092', None, u'ENSG00000233750')
# returns triple (upstream_gene, exact_gene, downstream_gene)
