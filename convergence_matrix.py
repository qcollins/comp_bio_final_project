from __future__ import print_function
import os
#import transcription_clusters
import yeast_emap_clusters
flatten = lambda l: [item for sublist in l for item in sublist]
indir = './100_out/'
#flat = flatten(flatten(transcription_clusters.gene_clusters))
flat = flatten(flatten(yeast_emap_clusters.gene_clusters))
relations = [[0 for x in range(len(flat))] for y in range(len(flat))] 
for root, dirs, filenames in os.walk(indir):
    for f in filenames:
        #clusters = file('out/'+f).readlines()
        clusters = file('100_out/'+f).readlines()
        for line in clusters:
            for wordrow in line.split():
                if("BPM" not in wordrow):
                    for wordcol in line.split():
                        if("BPM" not in wordcol):
                            relations[flat.index(wordrow)][flat.index(wordcol)] += 1

out = [[0 for x in range(101)] for y in range (2)]
for i in range(101):
    out[0][i] = i
for row in range(len(flat)):
    for col in range(len(flat)):
        #print(relations[row][col], end=" ")
        out[1][relations[row][col]] += 1
    #print("")
print(out)

