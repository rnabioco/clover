import tempfile
import argparse

class StandardizeAFA:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
        self.conversion = {
            'I': 'A', '3': 'T', '!': 'T', '$': 'T', 'N': 'T', '2': 'T', '#': 'G', ')': 'T',
            'Q': 'G', '1': 'T', 'P': 'T', 'V': 'T', '⊄': 'G', 'S': 'T', '{': 'T', 'B': 'C',
            'M': 'C', 'ʆ': 'G', 'J': 'T', '9': 'G', '.': 'A', 'ʭ': 'T', 'ƕ': 'T', '>': 'C', '*': 'A'
        }

    def convert_anticodon(self, anticodon):
        return f"{''.join(self.conversion.get(char, char) for char in anticodon)}-{anticodon}"

    def clean_fasta_header(self, intermediate_file):
        with open(self.input_file, 'r') as infile, open(intermediate_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    parts = line.split('|')
                    type_, amino_acid, anticodon, location = parts[1], parts[2], parts[3].replace('U', 'T'), parts[5].strip().lower()
                    if amino_acid == "Ini":
                        amino_acid = "iMet"
                    anticodon_converted = self.convert_anticodon(anticodon.strip())
                    prefix = 'mito' if 'mitochondrion' in location else 'nuc'
                    outfile.write(f'>{prefix}-{type_}-{amino_acid}-{anticodon_converted}\n')
                else:
                    outfile.write(line)

    def clean_and_modify_sequences(self, intermediate_file):
        unique_entries = set()
        with open(intermediate_file, 'r') as infile, open(self.output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    header = line.strip()
                    if len(header.split('-')) < 4:
                        continue  # Skip invalid headers
                    sequence = next(infile).strip()
                    modified_sequence = f'CCUAAGAGCAAGAAGAAGCCUGGN{sequence}GGCUUCUUCUUGCUCUUAGGAAAAAAAAAA'
                    if (header, modified_sequence) not in unique_entries:
                        unique_entries.add((header, modified_sequence))
                        outfile.write(header + '\n' + modified_sequence + '\n')

    def process(self):
        with tempfile.NamedTemporaryFile(mode='w+', delete=True) as temp_file:
            self.clean_fasta_header(temp_file.name)
            temp_file.seek(0)
            self.clean_and_modify_sequences(temp_file.name)

# Argument Parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standardize AFA files")
    parser.add_argument("input_file", type=str, help="Path to input AFA file")
    parser.add_argument("output_file", type=str, help="Path to output AFA file")
    args = parser.parse_args()

    processor = StandardizeAFA(args.input_file, args.output_file)
    processor.process()
    print(f"Standardized AFA saved to {args.output_file}")