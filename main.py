from graph import DeBrujinGraph
import gzip
import time
import random
import string


def read_fastq(path):
    with gzip.open(path, 'rt') as f:
        for line in f:
            sequence = f.readline().rstrip()
            if sequence[0] != '~':
                yield sequence

def read_fasta(path) -> [str]:
    with gzip.open(path, 'rt') as f:
        accession, description, seq = None, None, None
        for line in f:
            if line[0] == '>':
                # yield current record
                if accession is not None:
                    yield accession, description, seq

                # start a new record
                accession, description = line[1:].rstrip().split(maxsplit=1)
                seq = ''
            else:
                seq += line.rstrip()

def read_fasta_uncompressed(path) -> [str]:
    f= open(path, "r")

    accession, description, seq = None, None, None
    for line in f:
        if line[0] == '>':
            # yield current record
            if accession is not None:
                yield accession, description, seq

            # start a new record
            accession, description = line[1:].rstrip().split(maxsplit=1)
            seq = ''
        else:
            seq += line.rstrip()

def generate_kmers(sequences, k):
    kmer=[]
    for sequence in sequences:
        l = len(sequence)
        for i in range(l - k + 1):
            #yield (sequence[i:i + k])
            kmer.append(sequence[i:i + k])
    return kmer

def save_contigs(contigs):
    f = open("contigs.fa", "w+")

    random.seed(123)

    for i in range(len(contigs)):
        f.write(">" + ''.join(random.choices(string.ascii_uppercase, k=10)) + " Description, description" "\n")
        f.write(contigs[i] + "\n")
    f.write(">EndID End, End\n")

def find_occurences():
    f = open("occurences.bed", "w+")

    for i in read_fasta("GCF_000002985.6_WBcel235_rna.fna.gz"):
        for j in read_fasta_uncompressed("contigs.fa"):
            occurence_index = i[2].find(j[2])
            if(occurence_index != -1):
                f.write(i[0])
                f.write("\t")
                f.write(str(occurence_index+1))
                f.write("\t")
                f.write(str(occurence_index + 1 + len(j[2])))
                f.write("\t")
                f.write(j[0])
                f.write("\n")

def kmer_walk(kmer_graph, start):
    k = start
    # yield the starting node
    yield k
    while True:
        for symbol in ['A', 'T', 'C', 'G']:
            candidate = k[1:] + symbol
            if candidate in kmer_graph:
                k = candidate
                yield k
        else:
            break

def walk(graph,start,closed=set(),contig=None):
    for k in kmer_walk(graph, start=start):
        if contig is None:
            contig = k
        else:
            contig += k[-1]

        if k in closed:
            inner_walk(graph,k,closed)
            break
        else:
            closed.add(k)
    yield contig

# ce qui nous permet de continuer lorsque nous rencontrons un cycle
def inner_walk(graph,node,set,):
    succesors_list=[i for i in graph.succesors(node)]
    if succesors_list != []:
        for  successor in succesors_list:
            if successor not in set:
                walk(graph,node,set)
    else:
        pass

if __name__ == '__main__':
    start=time.time()
    #Question 2)a)

    sequences=read_fastq('reads.fastq.gz')
    graph=DeBrujinGraph(generate_kmers(sequences,21),21)

    #Question 3)a)
    starting_nodes=[i for i in graph.find_start()]
    #initialilisation du tableau fournissant les contigs
    contigs=[]

    #remplissage de contig dans le tableau pour chaque route emprunter
    for node in range (0,100):
        [contigs.append(i) for i in walk(graph,starting_nodes[node])]
    save_contigs(contigs)

    #Question 4)a)
    find_occurences()
    end=time.time()
    print(end-start, "secs")













