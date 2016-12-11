#from __future__ import print_function
import math
import random
from simanneal import Annealer
import gene_inter
from itertools import combinations, product
import csv
import prune

def normpdf(x, mean, sd):
    var = float(sd)**2
    pi = 3.1415926
    denom = (2*pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def interweight((A, B)):
    '''
    Calculates the interaction weight of a BPM. It is defined as the difference
    of sums of interaction scores within each module and the sum of interaction
    scores between each module divided by the number of genes in the entire BPM.

    The value returned is a BPM "decorated" with the interaction weight for
    sorting purposes. This roundabout means of decoration is used so that
    parallelism can be used for calculating the interaction weights. (As
    opposed to using a higher order function with 'sorted'.)
    '''

    if len(A) + len(B) == 0:
        return 0
    # For converting a tuple to two arguments
    gitup = lambda (g1, g2): gene_inter.gi(g1, g2)

    def sum_within(S):
        return sum(map(gitup, combinations(S, 2)))

    within = sum_within(A) + sum_within(B)
    between = sum(map(gitup, product(A, B)))


    iweight = (within - between) #/ float(len(A) + len(B))

    #print "A: ", A, "B ", B, "W: ", iweight, sum_within(A), sum_within(B), between
    return iweight

def scaled_interweight((A, B)):
    '''
    Calculates the interaction weight of a BPM. It is defined as the difference
    of sums of interaction scores within each module and the sum of interaction
    scores between each module divided by the number of genes in the entire BPM.

    The value returned is a BPM "decorated" with the interaction weight for
    sorting purposes. This roundabout means of decoration is used so that
    parallelism can be used for calculating the interaction weights. (As
    opposed to using a higher order function with 'sorted'.)
    '''
    gitup = lambda (g1, g2): gene_inter.gi(g1, g2)

    product_ab = filter(lambda x: gitup(x) != 0, list(product(A, B)))
    combo_a =    filter(lambda x: gitup(x) != 0, list(combinations(A, 2)))
    combo_b =    filter(lambda x: gitup(x) != 0, list(combinations(B, 2)))



    #print(len(product_ab), len(combo_a), len(combo_b))


    if len(A) + len(B) == 0:
        return 0
    # For converting a tuple to two arguments

    def sum_within(S):
        return sum(map(gitup, combinations(S, 2)))

    within = sum_within(A)/max(float(len(combo_a)), 1.0) + sum_within(B)/max(1.0, float(len(combo_b)))
    between = sum(map(gitup, product(A, B))) / max(float(len(product_ab)), 1.0)

    iweight = (within - between) #/ float(len(A) + len(B))

    return iweight




def bpm_energy((A, B)):
    #return scaled_interweight((A, B))
    #print(len(A) + len(B))
    w = interweight((A, B)) * normpdf(len(A) + len(B), 30, 5)
    return  -w #if w < 0 else w * normpdf(len(A) + len(B), 15, 5)

# State is (energy,[bpms]), where a bpm is 
#    ([nodes in module1], [nodes in module 2]) 

def randomBool():
    return bool(random.getrandbits(1))

def pop_random_gene(bpm):
    (module_1, module_2) = bpm
    if(randomBool() and len(module_1)!= 0 or not module_2):
        return module_1.pop(random.randrange(len(module_1)))
    elif(len(module_2) != 0):
        return module_2.pop(random.randrange(len(module_2)))
    else:
        print(bpm)
        raise NameError('poping gene from empty bpm')

def add_random_gene(gene, bpm):
    (module_1, module_2) = bpm
    if(randomBool()):
        module_1.append(gene)
    else:
        module_2.append(gene)

def is_bpm_empty((module_1, module_2)):
    return (not module_1 and not module_2)


class GeneAnnealer(Annealer):    
    # pass extra data (the distance matrix) into the constructor
    def __init__(self, state):
        #print(state)
        (energy,bpms) = state
        bpms = filter(lambda x: not is_bpm_empty(x), bpms)
        super(GeneAnnealer, self).__init__((energy, bpms))  # important! 
        

    def move(self):
        (energy, bpms) = self.state
        
        source_index = random.randrange(len(bpms)) 
        source_bpm = bpms[source_index]
        (a,b) = source_bpm
        if(len(a) > 25 or len(b) > 25):
            init_energy = bpm_energy(source_bpm)
            a_new = a[len(a)/2:]
            b_new = b[len(b)/2:]
            
            del a[len(a)/2:]
            del b[len(b)/2:]
            source_energy_new = bpm_energy(source_bpm)
            bpm_new = (a_new, b_new)
            bpm_new_energy = bpm_energy(bpm_new)
            bpms.append(bpm_new)
            energy = energy + bpm_new_energy + source_energy_new - init_energy
            self.state = (energy, bpms)
        else:
            target_index = random.randrange(len(bpms) + 1)
            if target_index == len(bpms):
                bpms.append(([], []))

            target_bpm = bpms[target_index]

            source_energy = bpm_energy(source_bpm)
            target_energy = bpm_energy(target_bpm)


            source_gene = pop_random_gene(source_bpm)
            add_random_gene(source_gene, target_bpm)
            
            source_energy_new = None
            if(is_bpm_empty(source_bpm)):
                del bpms[source_index]
                source_energy_new = 0
            else:
                source_energy_new = bpm_energy(source_bpm)
            target_energy_new = bpm_energy(target_bpm)
            energy_old = energy

            if target_index == source_index:
                energy = energy + target_energy_new - target_energy
            else:
                energy = (energy + target_energy_new + source_energy_new  
                             - target_energy - source_energy)
            
            # print "---------------------------------------"
            # print energy_old, energy, target_energy_new, source_energy_new, target_energy, source_energy
            # print "actual ", sum(map(bpm_energy, self.state[1])), "real ", energy, "state", self.state[1]
            self.state = (energy, bpms) 


    
    

    def energy(self):
        return self.state[0]
        #print sum(map(bpm_energy, self.state[1]))
        #return sum(map(bpm_energy, self.state[1]))


# TODO: loosing some genes in this process
def split_list(a_list):
    half = len(a_list)/2
    return a_list[:half], a_list[half:]

# partition to lists
def partition_set(set, num_partitions):
    partitions = [[] for i in range(num_partitions)]

    for elem in set:
        partitions[random.randrange(num_partitions)].append(elem)
    #print(partitions)
    return partitions

# partitions is list of lists
def pair_groups(partitions):
    (first, second) = split_list(partitions)
    return zip(first, second)


NUM_BPMS = 120


if __name__ == '__main__':
    out_fname = "out.bpm"#"results/norm_not_pruned.bpm" #"results/yeast_raw_unfiltered.bmps"
    gene_inter.load_genes("data/yeast_emap.gi", ignore_file="data/essentials")
    #gene_inter.load_genes("data/8_elem_test.gi")
    
    #print(gene_inter.gis)
    #print(gene_inter.genes)   

    initial_state = pair_groups(partition_set(gene_inter.genes, NUM_BPMS * 2)) 
    print(initial_state)
    energy = sum(map(bpm_energy, initial_state))
    annealer = GeneAnnealer((energy, initial_state))
    annealer.copy_strategy = "deepcopy"  
    state, e = annealer.anneal()

    print(state)

    (energy, bpms) = state

    # now prune bpms
    #bpms = prune.prune(bpms)

    print(bpms)
    outf = open(out_fname, 'w+')
    out = csv.writer(outf, delimiter='\t')
    for i, (mod1, mod2) in enumerate(bpms):
        mod1, mod2 = list(mod1), list(mod2)
        print("mod1: %d %d %d" %(len(mod1), len(mod2), bpm_energy((mod1, mod2))))
        out.writerow(['BPM%d/Module1' % i] + mod1)
        out.writerow(['BPM%d/Module2' % i] + mod2)
