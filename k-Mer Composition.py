from collections import defaultdict

def get_4mer_composition(dna_string):
    kmer_length = 4
    kmer_composition = defaultdict(int)

    # Use a sliding window to update the composition efficiently
    for i in range(len(dna_string) - kmer_length + 1):
        kmer = dna_string[i:i + kmer_length]
        kmer_composition[kmer] += 1

    # Output the counts in lexicographic order
    result = [kmer_composition.get(kmer, 0) for kmer in sorted(kmer_composition.keys())]

    return result

# Provided FASTA input
fasta_input = ">Rosalind_6431\nCTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG\nCCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT\nTTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCA\nAATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCG\nGGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGA\nCTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTA\nCCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG"

# Extract the DNA string from the FASTA input
dna_string = ''.join(line.strip() for line in fasta_input.split('\n')[1:])
result = get_4mer_composition(dna_string)
print(result)
