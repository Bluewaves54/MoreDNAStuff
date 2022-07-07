from random import choices

amino_acids = {
    'UUU': 'F', 'UUC': 'F',
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUA': 'L', 'CUG': 'L', 'CUC': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
    'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y',
    'CAU': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C',
    'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

stop_codons = ['UAA', 'UAG', 'UGA']

comps = {
    'A': 'U',
    'U': 'A',
    'C': 'G',
    'G': 'C'
}


def flatten(iters):
    return [elem for iter in iters for elem in iter]


def randseq(gcpercent=50, length=30000):
    sequence = ''.join(choices(
        population=['A', 'U', 'G', 'C'],
        weights=(100 - gcpercent, 100 - gcpercent, gcpercent, gcpercent),
        k=length
    ))
    return sequence


class Genome:
    def __init__(self, file_path=None, taxid=None, organism=None, genome=None, gcpercent=50, length=1000, name='unnamed'):
        if genome is not None:
            self.genome = genome
            self.name = name
        elif file_path is not None:
            with open(file_path) as file:
                lines = file.readlines()
                self.genome = ''.join([line.strip() for line in lines[1:]]).replace('T', 'U').replace('Y', 'U')
                self.name = lines[0][1:].replace('/', '-')
        else:
            self.genome = randseq(gcpercent=gcpercent, length=length)
        self.complement = ''.join([comps[let] for let in self.genome])
        self.rev_comp = self.complement[::-1]
        self.taxid = taxid
        self.organism = organism
        self.gc_content = (self.genome.count('G') + self.genome.count('C')) / self.genome.count(
            'A') + self.genome.count('U')

    # def __repr__(self):
    #     return self.genome

    def __getitem__(self, item):
        return self.genome[item]

    def get_protein_seqs(self, rfmode=3, start_codon=False, min_length=0):
        if rfmode == 3:
            protein_sequences = {i: [''] for i in range(3)}
            rflistnums = [0, 0, 0]
            for index, char in enumerate(self.genome[:-3:]):
                codon = self.genome[index:index + 3]
                rf = index % 3
                if start_codon:
                    if codon != 'AUG' and protein_sequences[rf][rflistnums[rf]] == '':
                        continue
                if codon in stop_codons:
                    rflistnums[rf] += 1
                    protein_sequences[rf].append('')
                else:
                    protein_sequences[rf][rflistnums[rf]] += amino_acids[codon]
            protein_sequences = {i: [j for j in protein_sequences[i] if len(j) >= min_length] for i in range(3)}
            protein_sequences = flatten(protein_sequences.values())
        if rfmode == 1:
            protein_sequences = ['']
            num = 0
            for index in range(0, len(self.genome[:-3:]), 3):
                codon = self.genome[index:index + 3]
                if start_codon:
                    if codon != 'AUG' and protein_sequences[num] == '':
                        continue
                if codon in stop_codons:
                    num += 1
                    protein_sequences.append('')
                else:
                    protein_sequences[num] += amino_acids[codon]
            protein_sequences = [j for j in protein_sequences if len(j) >= min_length]
        self.protein_seqs = protein_sequences
        return protein_sequences

    def get_proteins(self):
        """BLAST search protein sequences"""
        pass

    def get_organism(self):
        """get organism and taxid from arbitrary protein name"""
        pass
