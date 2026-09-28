import numpy
from sborf.predict import run_prediction
import numpy as np

import random

aa2cod = {'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'C': ['TGT', 'TGC'], 'E': ['GAG', 'GAA'], 'D': ['GAT', 'GAC'],
          'G': ['GGT', 'GGG', 'GGA', 'GGC'], 'F': ['TTT', 'TTC'], 'I': ['ATC', 'ATA', 'ATT'], 'H': ['CAT', 'CAC'],
          'K': ['AAG', 'AAA'], 'M': ['ATG'], 'L': ['CTC', 'CTG', 'CTA', 'CTT', 'TTA', 'TTG'], 'N': ['AAC', 'AAT'],
          'Q': ['CAA', 'CAG'], 'P': ['CCT', 'CCG', 'CCA', 'CCC'], 'S': ['AGC', 'AGT', 'TCG', 'TCA', 'TCC', 'TCT'],
          'R': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG', 'CGT'], 'T': ['ACC', 'ACA', 'ACG', 'ACT'], 'W': ['TGG'],
          'V': ['GTA', 'GTC', 'GTG', 'GTT'], 'Y': ['TAT', 'TAC']}

def make_combinations(prot_seq, ite=10):
    fs = []
    if ite > 0:
        for i in range(ite):
            s = []
            for k in prot_seq:
                s += random.sample(aa2cod[k], 1)
            fs += [s]
    return fs


cod2aa = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
          "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
          "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
          "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",

          "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
          "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
          "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
          "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",

          "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
          "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
          "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
          "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",

          "TAA": "_", "TAC": "Y", "TAG": "_", "TAT": "Y",
          "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
          "TGA": "_", "TGC": "C", "TGG": "W", "TGT": "C",
          "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}


def checkseq(seq, target):
    s = ''
    for i in seq:
        s += cod2aa[i]
    for i in range(len(s)):

        assert s[i] == target[i]
    assert len(s) == len(target)

    assert s == str(target)


def optimize(target_seq, organism, TARGET_SOL=1,iterations = 100,pop_size = 1500,num_parents_mating = 200,verbose=False,num_optimized_seqs=1,device="cpu",cancel_check=None):
    if verbose:
        print("\toptimization tith target sol ",TARGET_SOL)

    new_population = make_combinations(target_seq, ite=pop_size)

    for generation in range(iterations):

        if cancel_check is not None:
            cancel_check()  # may raise to abort the optimization

        fitness, preds = cal_pop_fitness(TARGET_SOL, new_population, organism,device=device)

        parents = select_mating_pool(new_population, fitness, num_parents_mating)

        offspring_crossover = crossover(parents, offspring_size=(pop_size - len(parents), len(parents[0])))

        offspring_mutation = mutation(offspring_crossover, target_seq)

        new_population = parents + offspring_mutation

        for i in new_population:
            checkseq(i, target_seq)

        best_match_idx = numpy.where(fitness == numpy.max(fitness))
        if verbose:
            print("gen:", generation, "fitness : ", fitness[best_match_idx[0][0]])

    fitness, pred = cal_pop_fitness(TARGET_SOL, new_population, organism,device=device)
    sorted_strings = ["".join(x) for _, x in sorted(zip(fitness, new_population),reverse=True)][:num_optimized_seqs]

    return sorted_strings


def cal_pop_fitness(target_sol, pop, organism,device):
    y = np.array(run_prediction(pop, organism,device=device))

    fitness =1 - (y - target_sol)**2
    return fitness, y


def select_mating_pool(pop, fitness, num_parents):
    return [x for _, x in sorted(zip(fitness, pop), reverse=True)][:num_parents]

def crossover(parents, offspring_size):
    offspring = []
    crossover_point = int(offspring_size[1] // 2)
    for k in range(offspring_size[0]):
        parent1_idx = k % len(parents)  # .shape[0]
        parent2_idx = (k + 1) % len(parents)

        offspring += [parents[parent1_idx][0:crossover_point] + parents[parent2_idx][crossover_point:]]

    return offspring


def mutation(offspring_crossover, target_seq, nmut=30):
    # Mutation changes a single gene in each offspring randomly.

    for idx in range(len(offspring_crossover)):
        r = np.random.randint(0, len(offspring_crossover[idx]), size=(nmut))
        for k in range(nmut):
            c = random.sample(aa2cod[target_seq[r[k]]], 1)[0]
            assert cod2aa[c] == target_seq[r[k]]

            offspring_crossover[idx][r[k]] = c
    return offspring_crossover



