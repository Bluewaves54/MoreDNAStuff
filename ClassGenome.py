from random import choices
from Bio.Seq import Seq

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
    """takes an iterable of iterables and returns a flattened iterable"""
    return [elem for iter in iters for elem in iter]


def randseq(gcpercent=50, length=30000):
    """
    creates a random RNA sequence of length 'length' (default is 30000),
    and a concentration of Gs and Cs of 'gcpercent' percent (default is 50%).
    """
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
                self.name = lines[0][1:].replace('/', '-').strip()
        else:
            self.genome = randseq(gcpercent=gcpercent, length=length)
            self.name = name
        self.complement = ''.join([comps[let] for let in self.genome])
        self.rev_comp = self.complement[::-1]
        self.taxid = taxid
        self.organism = organism
        self.gc_content = (self.genome.count('G') + self.genome.count('C')) / self.genome.count(
            'A') + self.genome.count('U')

    def __repr__(self):
        return self.name

    def __getitem__(self, item):
        return self.genome[item]

    def __eq__(self, other):
        if not isinstance(other, Genome):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.genome == other.genome

    def get_protein_seqs1(self, rfmode=3, start_codon=False, min_length=0, start_index=0):
        if rfmode == 3:
            protein_sequences = {i: [''] for i in range(3)}
            rflistnums = [0, 0, 0]
            for index, char in enumerate(self.genome[start_index:-3]):
                codon = self.genome[index + start_index:index + start_index + 3]
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
            for index in range(0, len(self.genome[start_index:-3:]), 3):
                codon = self.genome[index + start_index:index + start_index + 3]
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

    def get_protein_seqs2(self, rfmode=3, start_codon=False, min_length=0, start_index=0):
        seq = Seq(self.genome[start_index:])
        translations = []
        for i in range(rfmode):
            if start_codon:
                translations.append([i.split('M', 1)[1] for i in seq[i:].translate().split('*') if len(i.split('M', 1)[1]) >= min_length ])
            else:
                translations.append([i for i in seq[i:].translate().split('*') if len(i) >= min_length])
        translations = [str(i) for i in flatten(translations)]
        self.protein_seqs = translations
        return translations

    # def get_protein_seqs3(self, rfmode=3, start_codon=False, min_length=0, start_index=0):

    def get_genes(self, rfmode=3, start_codon=False, min_length=0, start_index=0):
        genes = []
        rfs_full_seqs = [self.get_translation(start_index=i) for i in range(rfmode)]
        rfs_indiv_seqs = [[j for j in frame.split('*') if len(j) > min_length] for frame in rfs_full_seqs]
        for i in range(rfmode):
            current_full_seq = rfs_full_seqs[i]
            current_indiv_seqs = rfs_indiv_seqs[i]
            current_rf_genes = []
            for j in current_indiv_seqs:
                protein_index = current_full_seq.find(j)
                genome_index = protein_index * 3
                gene_length = len(j) * 3
                gene = self.genome[i + genome_index:i + gene_length + genome_index]
                current_rf_genes.append((j, gene))
            genes.append(current_rf_genes)
        genes = flatten(genes)
        self.genes = genes
        return genes

    def get_translation(self, start_index=0):
        iterable = self.genome[start_index:-3]
        protein_sequence = ''
        for i in range(0, len(iterable), 3):
            codon = self.genome[start_index + i:start_index + i + 3]
            if codon not in stop_codons:
                protein_sequence += amino_acids[codon]
            else:
                protein_sequence += '*'
        self.translation = protein_sequence
        return self.translation

    def get_proteins(self):
        """BLAST search protein sequences"""
        pass

    def get_organism(self):
        """get organism and taxid from arbitrary protein name"""
        pass
