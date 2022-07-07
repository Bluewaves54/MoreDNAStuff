from random import choices
from ClassGenome import Genome
from matplotlib import pyplot as plt

def randseq(gcpercent=50, length=30000):
    sequence = ''.join(choices(
        population=['A', 'U', 'G', 'C'],
        weights=(100-gcpercent, 100-gcpercent, gcpercent, gcpercent),
        k=length
    ))
    return Genome(genome=sequence)


def testPercentage(bps=30000, gcpercent=50, loops=1000, min_length=100):
    counts = []
    for i in range(loops):
        sequence = randseq(gcpercent=gcpercent, length=bps)
        long_protein_seqs = sequence.get_protein_seqs(min_length=min_length)
        counts.append(len(long_protein_seqs))
    mean = sum(counts) / loops
    return mean

x = [i for i in range(101)]
y = [testPercentage(gcpercent=i, loops=5) for i in range(101)]

plt.plot(x, y)
plt.show()
