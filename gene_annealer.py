from __future__ import print_function
import math
import random
from simanneal import Annealer
import gene_inter
from itertools import combinations, product
import csv

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

    iweight = (within - between) / float(len(A) + len(B))

    return iweight

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
        # print(len(bpms))
        
        source_index = random.randrange(len(bpms)) 
        source_bpm = bpms[source_index]
        source_energy = interweight(source_bpm)
        source_gene = pop_random_gene(source_bpm)

        target_index = random.randrange(len(bpms) + 1)
        if target_index == len(bpms):
        	bpms.append(([], []))

        target_bpm = bpms[target_index]
        target_energy = interweight(target_bpm)

        add_random_gene(source_gene, target_bpm)
        
        if(is_bpm_empty(source_bpm)):
            del bpms[source_index]
        
        source_energy_new = interweight(source_bpm)
        target_energy_new = interweight(target_bpm)
        energy = (energy - target_energy_new - source_energy_new  
                         + target_energy + source_energy)
        
        self.state = (energy, bpms) 


    
    def move2(self):
        # remove a gene from a random bpm. 
        gene_to_move = None
        from_bpm = None

        while not gene_to_move:
            from_bpm = self.state[1][random.randint(0, len(self.state[1]) - 1)]
            energy_from_old = interweight(from_bpm)
            (from_module_1, from_module_2) = from_bpm
            if randomBool() and len(from_module_1) != 0:
                gene_to_move = from_module_1.pop(random.randrange(len(from_module_1)))
            elif len(from_module_2) != 0:
                gene_to_move = from_module_2.pop(random.randrange(len(from_module_2)))
            energy_from_new = interweight(from_bpm)

        # randomly move a gene from one bpm to another
        target_bpm   = self.state[1][random.randint(0, len(self.state[1]) - 1)]
        energy_target_old = interweight(target_bpm)
        (target_module_1, target_module_2) = target_bpm
        if randomBool():
            target_module_1.append(gene_to_move)
        else:
            target_module_2.append(gene_to_move)
        energy_target_new = interweight(target_bpm)

        (energy, bpms) = self.state
        energy = energy - energy_target_new - energy_from_new + energy_target_old + energy_from_old
        self.state = (energy, bpms) 

    def energy(self):
        return self.state[0]
        #return -sum(map(interweight, self.state))


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


NUM_BPMS = 60


if __name__ == '__main__':
    out_fname = "output.bmps"
    gene_inter.load_genes("data/yeast_emap.gi")
#    gene_inter.load_genes("data/8_elem_test.gi")
    
    #print(gene_inter.gis)
    #print(gene_inter.genes)   

    initial_state = pair_groups(partition_set(gene_inter.genes, NUM_BPMS * 2)) 
    print(initial_state)
    energy = -sum(map(interweight, initial_state))
    annealer = GeneAnnealer((energy, initial_state))
    annealer.copy_strategy = "deepcopy"  
    state, e = annealer.anneal()

    print(state)

    (energy, bpms) = state
    print(bpms)
    outf = open(out_fname, 'w+')
    out = csv.writer(outf, delimiter='\t')
    for i, (mod1, mod2) in enumerate(bpms):
        mod1, mod2 = list(mod1), list(mod2)
        out.writerow(['BPM%d/Module1' % i] + mod1)
        out.writerow(['BPM%d/Module2' % i] + mod2)
    # # initial state, a randomly-ordered itinerary
    # init_state = list(cities.keys())
    # random.shuffle(init_state)

    # # create a distance matrix
    # distance_matrix = {}
    # for ka, va in cities.items():
    #     distance_matrix[ka] = {}
    #     for kb, vb in cities.items():
    #         if kb == ka:
    #             distance_matrix[ka][kb] = 0.0
    #         else:
    #             distance_matrix[ka][kb] = distance(va, vb)

    # tsp = TravellingSalesmanProblem(init_state, distance_matrix)
    # # since our state is just a list, slice is the fastest way to copy

    # while state[0] != 'New York City':
    #     state = state[1:] + state[:1]  # rotate NYC to start
    # print("%i mile route:" % e)
    # for city in state:
    #     print("\t", city)

