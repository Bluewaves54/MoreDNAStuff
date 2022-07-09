from Bio import SeqIO

def findInsertion():
    fasta_sequences = SeqIO.parse(open('spikeProteins.afa'), 'fasta')
    proteins = [(seq.id, str(seq.seq)) for seq in fasta_sequences]
    print(proteins)

    sarscov2_id = 'NC_045512'
    sars_protein = [i[1] for i in proteins if i[0] == sarscov2_id][0]
    print(sars_protein)

    inserted_codons = []
    inserted_indices = []
    suspected_genomes = ['MN994467', 'MN996532.2', 'MZ571142.1', 'Pangolin']

    current_frames = []
    for sars_index, char in enumerate(sars_protein[:-4:]):
        sars_frame = sars_protein[sars_index:sars_index + 4]
        for (index, (id, sequence)) in enumerate(proteins):
            if id != sarscov2_id:
                current_frames.append((id, sequence[sars_index:sars_index + 4]))

        set_current_frames = list(set(current_frames))
        sars_equal_frames = [frame[0] for frame in set_current_frames if frame[1] == sars_frame]
        print(sars_equal_frames)

        if set(sars_equal_frames) == set(suspected_genomes):
            inserted_codons.append(sars_frame)
            inserted_indices.append(sars_index)
        else:
            current_frames = []

    return inserted_indices, inserted_codons


def findCodonFrequencies():
    fasta_sequences = SeqIO.parse(open('spikeProteinSeqs.fasta'), 'fasta')
    proteins = ''.join([str(seq.seq) for seq in fasta_sequences])
    letters = {letter: proteins.count(letter) / len(proteins) * 100 for letter in proteins}
    four_letter_words = [proteins[i:i + 4] for i in range(len(proteins[:-4:]))]
    four_letter_word_frequencies = {word: proteins.count(word) / len(proteins) * 100 for word in four_letter_words}
    print(four_letter_word_frequencies)
    print(proteins.count('PRRA') / len(four_letter_words) * 100)

findCodonFrequencies()
