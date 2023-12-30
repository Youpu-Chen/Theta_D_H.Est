"""
Description: 
- fill the haplotype for the haploid-level theta calculation in Theta_D_F_H
- filter "." loci

Notes: Theta_D_F_H only allow diploid calculation by default, this script fills the haplotype and you can just specify the haplotype index (e.g., 1 or 2) for the calculation
Output: .vcf by default
"""
import gzip
import sys

def process_vcf_gz(input_file, output_file):
    with gzip.open(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # skip header lines
            if line.startswith('#'):
                outfile.write(line)
                continue

            # split the line into columns
            columns = line.strip().split('\t')

            # skip the line if it contains '.'
            if '.' in columns:
                continue

            # Transform the allele data
            for i in range(9, len(columns)):
                allele_data = columns[i].split(':')
                if allele_data[0] == '0':
                    allele_data[0] = '0|0'
                elif allele_data[0] == '1':
                    allele_data[0] = '1|1'
                columns[i] = ':'.join(allele_data)

            # Write the transformed line to the output file
            outfile.write('\t'.join(columns) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py input.vcf.gz output.vcf")
        sys.exit(1)

    input_vcf_gz = sys.argv[1]
    output_vcf = sys.argv[2]
    process_vcf_gz(input_vcf_gz, output_vcf)
