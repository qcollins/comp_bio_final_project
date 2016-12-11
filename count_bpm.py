from bpm import bpmreader
import sys
#
# Usage. Return
#

if len(sys.argv) < 2:
    print("Usage")
    exit()

file_name = sys.argv[1]
bpms = bpmreader.read(file_name)

flat = reduce(lambda rest, (a, b): a + b + rest, bpms, [])

print "Number genes: " + str(len(flat))

