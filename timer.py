import time
from ClassGenome import Genome
from Bio.Seq import Seq

genome = Genome(file_path='SARS2.fasta')

sequence = Seq(genome.genome)

start1 = time.perf_counter()
for i in range(500):
    protein_sequences = genome.get_protein_seqs1(min_length=50)
print(f'Genome class finished in {time.perf_counter() - start1}')

start2 = time.perf_counter()
for i in range(500):
    protein_sequences = genome.get_protein_seqs2(min_length=50)
print(f'Seq class finished in {time.perf_counter() - start2}')
