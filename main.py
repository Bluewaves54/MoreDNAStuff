import time
from CompareGenomes import compareGenomes
from ClassGenome import Genome, flatten, amino_acids, stop_codons
from ClassGenome import flatten
from os import listdir
from collections import Counter
from matplotlib import pyplot as plt


def getGenomesFromeDir(dirpath):
    fasta_files = [file for file in listdir(dirpath) if file.endswith('fasta')]
    genomes = [Genome(file_path=f'{dirpath}/{i}') for i in fasta_files]

def alignProteinGene(protein, gene):
    print(*list(protein), sep='\t')
    print(*[gene[i:i + 3] for i in range(0, len(gene), 3)])

def separateGenomeIntoCodons(gene):
    return flatten([[gene[i + j:i + j + 3] for i in range(0, len(gene), 3)] for j in range(3)])


if __name__ == "__main__":
    # start = time.perf_counter()
    # print(*compareGenomes('TEST_DIR'), sep='\n')
    # print(f'Finished in {time.perf_counter() - start} seconds')
    sars_spike_protein = 'LEKTTELLFLVMFLLTTKRTMFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT'
    print(sars_spike_protein.find('PRRA'))
    # genome = Genome(file_path='SARS-class/SARS-2020.fasta')
    # genome.get_genes(rfmode=3, min_length=100)
    # special_codons = 'CCU CGG CGG GCA'.split(' ')
    # cgg = special_codons[1]
    # sars_spike_protein_gene = [i for i in genome.genes if i[0] == sars_spike_protein][0]
    # alignProteinGene(*sars_spike_protein_gene)
    # codons = separateGenomeIntoCodons(genome.genome)
    # real_amino_acids = [i for i in codons if i not in stop_codons and len(i) == 3]
    # frequencies = {codon: frequency / len(real_amino_acids) * 100 for codon, frequency in dict(Counter(real_amino_acids)).items()}
    # print(frequencies)
    # r_frequencies = {key: value for key, value in frequencies.items() if amino_acids[key] == 'R'}
    # print(dict(sorted(r_frequencies.items(), key=lambda item: item[1])))
    # print(separateGenomeIntoCodons(genome.genome))


