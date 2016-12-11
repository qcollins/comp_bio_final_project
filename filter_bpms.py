from bpm import bpmreader
import sys
import csv
#
# Usage. Return
#
MAX = 25
MIN = 3


if len(sys.argv) < 3:
    print("Usage")
    exit()


if len(sys.argv) == 5:
    MIN = int(sys.argv[3])
    MAX = int(sys.argv[4])

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]

bpms = bpmreader.read(input_file_name)

def filter_fun((A, B)):
    a_len = len(A)
    b_len = len(B)

    return (a_len <= MAX and a_len >= MIN and
            b_len <= MAX and b_len >= MIN)



bpms = filter(filter_fun, bpms)


outf = open(output_file_name, 'w+') if output_file_name != "-" else sys.stdout
out = csv.writer(outf, delimiter='\t')
for i, (mod1, mod2) in enumerate(bpms):
    mod1, mod2 = list(mod1), list(mod2)
    out.writerow(['BPM%d/Module1' % i] + mod1)
    out.writerow(['BPM%d/Module2' % i] + mod2)
