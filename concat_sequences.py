import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus

def read_sequences(input_files, file_format):
    sequences = {}
    for input_file in input_files:
        if file_format == 'nexus':
            with open(input_file) as handle:
                nexus = Nexus.Nexus(handle)
                for taxon, seq in nexus.data.items():
                    if taxon not in sequences:
                        sequences[taxon] = []
                    sequences[taxon].append(seq)
        else:
            for record in SeqIO.parse(input_file, file_format):
                if record.id not in sequences:
                    sequences[record.id] = []
                sequences[record.id].append(str(record.seq))
    return sequences

def concatenate_sequences(sequences, fill_missing, sequence_length):
    all_taxa = set(sequences.keys())
    concatenated = {}

    for taxon, seqs in sequences.items():
        concatenated[taxon] = ''.join(seqs)

    if fill_missing:
        max_length = max(len(seq) for seq in concatenated.values())
        for taxon in all_taxa:
            if taxon not in concatenated:
                concatenated[taxon] = 'N' * sequence_length

    return concatenated

def write_output(concatenated, output_file):
    with open(output_file, 'w') as output_handle:
        for taxon, sequence in concatenated.items():
            output_handle.write(f">{taxon}\n{sequence}\n")

def main():
    parser = argparse.ArgumentParser(description="Concatenate sequence data from FASTA, Phylip, or Nexus formats.")
    parser.add_argument('input_files', nargs='+', help='Input sequence files')
    parser.add_argument('-f', '--format', required=True, choices=['fasta', 'phylip', 'nexus'], help='Input file format')
    parser.add_argument('-o', '--output', required=True, help='Output concatenated file')
    parser.add_argument('--fill-missing', action='store_true', help='Fill missing sequences with Ns')
    parser.add_argument('--sequence-length', type=int, default=0, help='Length of the sequences to fill with Ns (if fill-missing is used)')
    
    args = parser.parse_args()

    sequences = read_sequences(args.input_files, args.format)
    concatenated = concatenate_sequences(sequences, args.fill_missing, args.sequence_length)
    write_output(concatenated, args.output)

if __name__ == '__main__':
    main()
