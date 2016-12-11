from bpm import enrichment
import sys
import gene_inter
import gene_annealer
from collections import defaultdict, OrderedDict


if len(sys.argv) < 2:
    print("Usage")
    exit()

input_file_name = sys.argv[1]

bpms = defaultdict(dict)

for bpmtext in open(input_file_name).read().split('>')[1:]:
        bpmi, modi, genes, goterms = enrichment.read_bpm(bpmtext)
        #print  goterms.keys()
        bpms[bpmi][modi] = (goterms.keys(), genes)
        #enriched[bpmi][modi] = None

gene_inter.load_genes("data/yeast_emap.gi", ignore_file="data/essentials")


for v in bpms.values():

    (go0, mod0) = v[0]
    (go1, mod1) = v[1]



    print "(%d  = %d & %d)\t%d\t%d\t%f" % (len(set(go0) & set(go1)), len(go0), len(go1), len(mod0), len(mod1), gene_annealer.bpm_energy((mod0, mod1)))
# print bpms
