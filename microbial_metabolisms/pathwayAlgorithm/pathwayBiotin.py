## This program takes reformatted KEGG Definitions for a pathway 
## (from the links below) and converts it into a graph.
#  
## I take the start and end KEGG entries and use them as source and 
## target nodes respectively. 
# 
## I generate all possible paths from source and target nodes and 
## check whether a MAG contains all the genes in a pathway with up to
## THRESHOLD number of genes missing. 
#
## THRESHOLD is my margin of error to account for genes KEGG may not have 
## picked up

## KEGG definitions from:
## https://www.genome.jp/module/M00123



from platform import node
from igraph import *

THRESHOLD = 2

gene_to_node = {'K00652':1,'K00833':2,'K19563':3,'K19562':4,'K01935':5,'K01012':6,'K16593':7,'K01906':8}
node_to_gene = dict((v,k) for k,v in gene_to_node.items())
source_vertices = ['K00652','K16593','K01906']
target_vertices = ['K01012',]
edges = [
    (1,2), (1,3), (1,4), 
    (2,5), (2,5),
    (5,6), (4,6), 
    (7,1), (1,3), (3,5), (5,6),
    (8,1), (1,2), (1,3), 
    (2,5), (3,5),
    (5,6)
]
g = Graph(edges, directed=True)

#setting up graph for visualizing
for i in list(node_to_gene.keys()):
    g.vs[i]["label"] = node_to_gene[i]
# plot(g)

#setting up list of graphs 
listOfPathways = list()
for s in source_vertices:
    for t in target_vertices:
        a = g.get_all_simple_paths(gene_to_node[s],gene_to_node[t])
        for pathway in a:
            listOfPathways.append(pathway)
print(listOfPathways)

def checkFx(magSet):
    output = []
    for pathway in listOfPathways:
        missing = 0
        count = 0
        ignore = False
        for node in pathway:
            count += 1
            if node_to_gene[node] in magSet:
                pass
            else:
                # print(node_to_gene[node])
                missing += 1
                if missing == THRESHOLD:
                    ignore = True
                    break
        
        if (ignore == False) and (len(pathway) > 1):
            output.append((missing, count, [node_to_gene[node] for node in pathway]))
    
    if output == []:
        return False
    else:
        return True

basepath = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/masterdata_17May2021"

# Load MAG names
misc_list = {}
with open(f"{basepath}/geneFunctionTables/goodBins/tatoosh_MAG_names.txt",'r') as m:
    m = m.readlines()
    for mag in m:
        misc_list[mag.strip()] = []

# Set up masterdata sheet
mag_gene_dict = {}
with open(f"{basepath}/geneFunctionTables/goodBins/tatoosh_masterdata.txt", 'r') as g:
    g = g.readlines()
    for line in g:
        line = line.split()
        if float(line[-2]) < 1e-20: # COG hits under 1e-50
            if line[3] == "KOfam":
                if line[2] in mag_gene_dict:
                    if "!!!" in line[4]:
                        for cog in line[4].split("!!!"):
                            mag_gene_dict[line[2]].append(cog)
                    else:
                        mag_gene_dict[line[2]].append(line[4])
                else: # First time mag name appears
                    if "!!!" in line[4]:
                        mag_gene_dict[line[2]] = line[4].split("!!!")
                    else:
                        mag_gene_dict[line[2]] = [line[4]]

output = ["\t".join(["MAG_names","env","Presence"])]

for mag in misc_list:
    sampleID = "_".join(mag.split("_")[:2])

    if checkFx(mag_gene_dict[mag])==True:
        output.append("\t".join([mag, sampleID, "Y"]))
    else:
        output.append("\t".join([mag, sampleID, "N"]))

with open(f"{basepath}/heatmapGen/heatmapFiles/pathwayBiotin.txt", 'w') as n:
    n.write("\n".join(output))




