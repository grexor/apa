import apa
import os
import sys
import pybio

def count_table(exp_list, output=""):

    species = set()
    expression = {}
    names = []
    for (lib_id, exp_id) in exp_list:
        lib = apa.annotation.libs[lib_id]
        exp_data = lib.experiments[exp_id]
        species.add(exp_data["map_to"])
        k = "%s_%s [%s]" % (lib_id, exp_id, exp_data["condition"])
        e_filename = apa.path.e_filename(lib_id, exp_id)
        expression[k] = pybio.data.Bedgraph(e_filename)
        names.append(k)

    assert(len(species)==1)
    species = list(species)[0]

    # find common list of positions
    positions = {}
    for id, bg in expression.items():
        for chr, strand_data in bg.raw.items():
            positions.setdefault(chr, {})
            for strand, pos in strand_data.items():
                positions.setdefault(chr, {}).setdefault(strand, set())
                positions[chr][strand] = positions[chr][strand].union(set(pos.keys()))

    # organize polya sites inside genes
    genes_sites = {}
    for chr, strand_data in positions.items():
        for strand, pos_set in strand_data.items():
            pos_set = list(pos_set)
            for pos in pos_set:
                gid_up, gid, gid_down, gid_interval = apa.polya.annotate_position(species, chr, strand, pos)
                if gid==None:
                    continue# only consider polya sites inside genes
                sites = genes_sites.get(gid, [])
                cDNA_sum = 0
                expression_vector = []
                for exp_id in names:
                    bg = expression[exp_id]
                    cDNA = bg.get_value(chr, strand, pos)
                    cDNA_sum += cDNA
                    expression_vector.append(cDNA)
                site_record = {"cDNA_sum":int(cDNA_sum), "site_pos":int(pos)} # later list of dictionaries will be sorted, still works fine (cDNA is first position, pos is second position)
                sites.append(site_record)
                genes_sites[gid] = sites

    # expression.genes
    genes_filename = output+"genes.tab"
    f_genes = open(genes_filename, "wt")
    header = ["chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype"]
    header += ["sites_position_cDNA", "sites_num", "cDNA_sum"]
    for exp_id in names:
        header.append(exp_id)
    f_genes.write("\t".join(header)+"\n")
    for gid, sites in genes_sites.items():
        sites = [(site_data["site_pos"], site_data["cDNA_sum"]) for site_data in sites]
        sites.sort()
        gene = apa.polya.get_gene(species, gid)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)
        row = [chr, strand, gene_locus, gid, gene["gene_name"], gene["gene_biotype"]]
        row.append(sites)
        row.append(len(sites))
        row_comp = []
        cDNA_total = 0
        for exp_id in names:
            bg = expression[exp_id]
            cDNA_comp = 0
            for (site_pos, _) in sites:
                val = bg.get_value(chr, strand, site_pos)
                cDNA_comp += val
            cDNA_total += cDNA_comp
            row_comp.append(cDNA_comp)
        row.append(cDNA_total)
        row += row_comp
        f_genes.write("\t".join([str(x) for x in row]) + "\n")
    f_genes.close()

    # expression.sites
    fname = output + "sites.tab"
    f = open(fname, "wt")
    header = ["site_id", "chr", "strand", "gene_locus", "gene_id", "gene_name", "gene_biotype", "site_pos"]
    for exp_id in names:
        header.append(exp_id)
    f.write("\t".join(header)+"\n")
    for gid, sites in genes_sites.items():
        sites = [(site_data["site_pos"], site_data["cDNA_sum"]) for site_data in sites]
        sites.sort()
        gene = apa.polya.get_gene(species, gid)
        chr = gene["gene_chr"]
        strand = gene["gene_strand"]
        gene_start = gene["gene_start"]
        gene_stop = gene["gene_stop"]
        gene_locus = "chr%s:%s-%s" % (chr, gene_start, gene_stop)
        for (site_pos, _) in sites:
            row = ["%s_%s" % (gid, site_pos), chr, strand, gene_locus, gid, gene["gene_name"], gene["gene_biotype"], site_pos]
            for exp_id in names:
                bg = expression[exp_id]
                val = bg.get_value(chr, strand, site_pos)
                row.append(val)
            f.write("\t".join([str(x) for x in row]) + "\n")
    f.close()
