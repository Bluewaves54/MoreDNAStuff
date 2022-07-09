from Bio.pairwise2 import format_alignment, align
from ClassGenome import Genome
from os import listdir
from itertools import combinations
from multiprocessing import Process, Manager
import time
import numpy as np

def compareGenomesChunks(genome_pairings, genome_score_avgs):
    scores = []
    for genome1, genome2 in genome_pairings:
        # print(genome1, genome2)
        for i, seq in enumerate(genome1.protein_seqs):
            for j, seq2 in enumerate(genome2.protein_seqs[i::]):
                alignment = align.globalxx(seq, seq2)
                scores.append(alignment)
        top_scores = []

        for i in range(len(genome1.protein_seqs)):
            top_scores.append(max(scores, key=lambda x: x[0][2] / len(x[0][1])))
            scores.remove(max(scores, key=lambda x: x[0][2] / len(x[0][1])))

        avg_score = sum([i[0][2] / len(i[0][1]) for i in top_scores]) / len(top_scores)

        with open(f'alignments/{genome1.name}x{genome2.name}.txt', 'a') as file:
            file.writelines([format_alignment(*i[0]) for i in top_scores])
        genome_score_avgs.append((genome1, genome2, avg_score))

def divideGenomes(dirpath):
    fasta_files = [file for file in listdir(dirpath) if file.endswith('fasta')]
    genomes = [Genome(file_path=f'{dirpath}/{i}') for i in fasta_files]
    for i in genomes:
        i.get_protein_seqs1(min_length=100, start_codon=False)

    genome_pairings = list(combinations(genomes, 2))
    genome_pairings_chunks = np.array_split(genome_pairings, 5)
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    return genome_pairings_chunks

def createTimeline(score_avgs):
    min_score = min(score_avgs, key=lambda x: x[2])
    min_pairs = [(i if i[0].genome == min_score[0].genome else (i[1], i[0], i[2])) for i in score_avgs]
    min_pairs = [i for i in score_avgs if min_score[0] == i[0]]
    min_pairs_sorted = sorted(min_pairs, key=lambda x: x[2], reverse=True)
    flattened_final = [i for j in min_pairs_sorted for i in j[:2:]]
    unique_genomes = []
    final = []
    for i in flattened_final:
        if i.genome not in unique_genomes:
            unique_genomes.append(i.genome)
            final.append(i)
    return final

def compareGenomes(dir_path):
    with Manager() as manager:
        genome_score_avgs = manager.list()
        processes = [Process(target=compareGenomesChunks, args=(chunk, genome_score_avgs,)) for chunk in divideGenomes('TEST_DIR')]
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        return createTimeline(genome_score_avgs)
