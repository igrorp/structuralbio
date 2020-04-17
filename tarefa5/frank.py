
def read_fasta(filename):
    with open(filename) as file:
        return [(part[0], part[2].replace('\n', ''))
            for part in [entry.partition('\n')
                for entry in file.read().split('>')[1:]]]


dict = {}

for header, seq in read_fasta('zea_prots.fasta'):

    gene = header.split(' ')[3][5:]

    if gene in dict:

        if len(dict[gene]) < len(seq):

            dict[gene] = seq

    else:

        dict[gene] = seq


with open('out.fasta', 'w') as nf:

    for gene in dict:

        nf.write(f'>{gene}\n{dict[gene]}\n')