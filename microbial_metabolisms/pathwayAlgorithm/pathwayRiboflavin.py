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
## https://www.genome.jp/module/M00911
## https://www.genome.jp/module/M00125

from platform import node
from igraph import *

THRESHOLD = 1

gene_to_node = {
    'K01497':1,'K14652':2,'K01498':3,'K00082':4,'K11752':5,
    'K22912':6,'K20860':7,'K20861':8,'K20862':9,'K21063':10,
    'K21064':11,'K20884':12,'K22949':13,'K11753':14,'K14654':15,
    'K14655':16,'K02858':17,'K00794':18,'K00793':19,'K00861':20,'K00953':21}
node_to_gene = dict((v,k) for k,v in gene_to_node.items())
source_vertices = ['K01497','K14652','K02858','K00794']
target_vertices = ['K22912','K20860','K20861','K20862','K21063','K21064', 'K02858','K14652', 'K22949','K11753', 'K14655','K00953']
g = Graph([(1,3),(1,5),(2,3),(2,5),(3,4),
(4,6),(4,7),(4,8),(4,9),(4,10),(4,11),
(5,6),(5,7),(5,8),(5,9),(5,10),(5,11),
(18,19),(19,12),(19,14),(12,13),
(1,15),(15,16),(19,20),(20,21)], directed=True)

#setting up graph for visualizing
for i in list(node_to_gene.keys()):
    g.vs[i]["label"] = node_to_gene[i]

plot(g)

#setting up list of graphs 
listOfPathways = list()
for s in source_vertices:
    for t in target_vertices:
        a = g.get_all_simple_paths(gene_to_node[s],gene_to_node[t])
        for pathway in a:
            listOfPathways.append(pathway)

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

with open(f"{basepath}/heatmapGen/heatmapFiles/pathwayRiboflavin.txt", 'w') as n:
    n.write("\n".join(output))




