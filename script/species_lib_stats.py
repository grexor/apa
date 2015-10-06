import apa

reads = {}

for lib_id, lib_data in apa.annotation.libs.items():
    fname = "data.apa/%s/%s_stats.tab" % (lib_id, lib_id)
    stats = {}
    f = open(fname, "rt")
    r = f.readline()
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        stats[int(r[0])] = float(r[-3])
        r = f.readline()
    f.close()
    for exp_id, exp_data in lib_data.experiments.items():
        if stats.get(exp_id, None)!=None:
            reads[exp_data["map_to"]] = reads.get(exp_data["map_to"], 0) + stats[exp_id]

for species, readM in reads.items():
    print "%s\t%s" % (species, readM)
