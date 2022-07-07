from Bio.pairwise2 import format_alignment, align
from ClassGenome import Genome
from os import listdir
from itertools import combinations

def compareMultipleGenomes(dirpath):
    fasta_files = [file for file in listdir(dirpath) if file.endswith('fasta')]
    print(fasta_files)
    genomes = [Genome(file_path=f'{dirpath}/{i}') for i in fasta_files]
    for i in genomes:
        i.get_protein_seqs(min_length=100, start_codon=False)

    genome_pairings = list(combinations(genomes, 2))

    genome_score_avgs = []
    for genome1, genome2 in genome_pairings:
        print(genome1.name, genome2.name)
        scores = []
        for i, seq in enumerate(genome1.protein_seqs):
            for j, seq2 in enumerate(genome2.protein_seqs[i::]):
                alignment = align.globalxx(seq, seq2)
                scores.append(alignment)
        top_scores = []

        for i in range(len(genome1.protein_seqs)):
            top_scores.append(max(scores, key=lambda x: x[0][2] / len(x[0][1])))
            scores.remove(max(scores, key=lambda x: x[0][2] / len(x[0][1])))

        avg_score = sum([i[0][2] / len(i[0][1]) for i in top_scores]) / len(top_scores)
        genome_score_avgs.append((genome1, genome2, avg_score))

        with open(f'alignments/{genome1.name}x{genome2.name}.txt', 'a') as file:
            file.writelines([format_alignment(*i[0]) for i in top_scores])
        print(genome_score_avgs)


compareMultipleGenomes('SARS-class')
