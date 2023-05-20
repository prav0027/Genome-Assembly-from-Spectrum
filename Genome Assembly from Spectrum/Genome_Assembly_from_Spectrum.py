from itertools import filterfalse
import networkx as nx

def read_reads(read_fn):
    all_reads = {}
    with open(read_fn, 'r') as f:
        #next(f)
        for line in f:
            if '>' in line:
                #line = line.strip()
                #header = line
                header = line[0:len(line)-1]
                all_reads[header] = ''
                continue
            #line = line.strip()
            #all_reads.append(line)
            else:
                #line = line.strip()
                #all_reads[header] = line
                all_reads[header] = line[0:len(line)-1]
    return all_reads

def debruijn_from_reads(reads, k):
    G = nx.MultiDiGraph()
    for read_id, read in reads.items():
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            G.add_edge(kmer[:-1], kmer[1:], readid=read_id)            
    return G

def save_results(sortedReads, fileName = 'predictions.csv'):
    f = open(fileName, 'w')
    for i in sortedReads:
        f.write(i+'\n')
    return

reads = read_reads('project3a_10000_spectrum.fasta')
graph = debruijn_from_reads(reads, 20)
path = list(nx.eulerian_path(graph))

# And map back to read headers
mapping = nx.get_edge_attributes(graph, "readid")
mappingTemp = dict()
for iKey, iValue in mapping.items():
    tupleTemp = (iKey[0], iKey[1])
    mappingTemp[tupleTemp] = iValue

print(mappingTemp)

sorted_reads = [mappingTemp[edge] for edge in path]
sortedReadsList = list()

for edge in path:
    print(edge)
    sortedReadsList.append(mappingTemp[edge])

save_results(sortedReadsList)

