from bpm import bpmreader
import sys
import csv


def filter_fun(MIN, MAX):
    def filterer((A, B)):

        a_len = len(A)
        b_len = len(B)

        return (a_len <= MAX and a_len >= MIN and
                b_len <= MAX and b_len >= MIN)
    return filterer




def filter_file(from_name, to_name, MIN, MAX):


    bpms = bpmreader.read(from_name)
    new_bpms = filter(filter_fun(MIN, MAX), bpms)



    outf = open(to_name, 'w+') if to_name != "-" else sys.stdout
    out = csv.writer(outf, delimiter='\t')
    for i, (mod1, mod2) in enumerate(new_bpms):
        mod1, mod2 = list(mod1), list(mod2)
        out.writerow(['BPM%d/Module1' % i] + mod1)
        out.writerow(['BPM%d/Module2' % i] + mod2)




if __name__ == '__main__':
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

    filter_file(input_file_name, output_file_name, MIN, MAX)
