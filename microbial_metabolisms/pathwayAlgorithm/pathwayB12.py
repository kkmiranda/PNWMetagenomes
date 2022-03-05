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
## https://www.genome.jp/module/M00924
## https://www.genome.jp/module/M00122
## https://www.genome.jp/module/M00925


from platform import node
from igraph import *


# Anaerobic
gene_to_node = {
    'K02303':1,'K02304':2,'K13542':3,'K24866':4,'K02302':5,
    'K02190':6,'K03795':7,'K22011':8,'K03394':9,'K05934':10,
    'K13541':11,'K21479':12,'K05936':13,'K02189':14,'K02188':15,
    'K05895':16,'K02191':17,'K00595':18,'K03399':19,'K06042':20,
    'K02224':21,'K13540':22,'K02229':23,'K02228':24,'K02230':25,
    'K09882':26,'K09883':27,'K00798':28,'K02232':29,'K19221':30,
    'K02225':31,'K02227':32,'K02231':33,'K00768':34,'K02226':35,'K22316':36
}
node_to_gene = dict((v,k) for k,v in gene_to_node.items())
source_vertices = ['K02303','K13542','K02302']
target_vertices = ['K02224','K09883']

anaerobicEdges = [
    #M00924 Anaerobic to cobyrinate a,c-diamide
    (1,2), (3,4), (1,4), (3,2),
    (5,6), (5,7), (5,8),
    (2,6), (2,7), (2,8),
    (4,6), (4,7), (4,8),
    (6,9), (7,9), (8,9),
    (9,10), (9,11), (9,12),
    (10,13), (11,13), (12,13),
    (13,14), (13,11),
    (14,15), (11,15),
    (15,16),
    (16,17), (16,18),
    (17,19),
    (19,20),(18,20),
    (20,21)]

aerobicEdges = [
    #M00925 Aerobic to cobyrinate a,c-diamide
    (1,9), (1,22),
    (3,9), (3,22),
    (9,23), (22,23),
    (23,10), (23,22), (23,11), 
    (10,13), (11,13),
    (13,24),
    (24,16),
    (16,18),
    (18,20),
    (20,21),
    (21,25),
    (25,26),
    (26,27)]

finalEdges = [
    # M00122 cobyrinate a,c-diamide => cobalamin
    (28,29), (30,29),
    (29,31), (29,32),
    (31,33), (32,33),
    (34,35), (34,36)
]

g_Anaerobic = Graph(anaerobicEdges, directed=True)
g_Aerobic = Graph(aerobicEdges, directed=True)
g_final = Graph(finalEdges, directed=True)

# #setting up graph for visualizing
# for i in list(node_to_gene.keys()):
#     g_Anaerobic.vs[i]["label"] = node_to_gene[i]
#     g_Aerobic.vs[i]["label"] = node_to_gene[i]
#     g_final.vs[i]["label"] = node_to_gene[i]


#setting up list of graphs 
anaerobicPathway = list()
for s in ['K02302', 'K02303', 'K13542']:
    for t in ['K02224']:
        a = g_Anaerobic.get_all_simple_paths(gene_to_node[s],gene_to_node[t])
        for pathway in a:
            anaerobicPathway.append(pathway)
aerobicPathway = list()
for s in ['K02303', 'K13542']:
    for t in ['K09883']:
        a = g_Aerobic.get_all_simple_paths(gene_to_node[s],gene_to_node[t])
        for pathway in a:
            aerobicPathway.append(pathway)

finalPathway = list()
for s in ['K00798', 'K19221', 'K00768']:
    for t in ['K02231', 'K02226', 'K22316']:
        a = g_final.get_all_simple_paths(gene_to_node[s],gene_to_node[t])
        for pathway in a:
            finalPathway.append(pathway)
pathwaySet = [anaerobicPathway, aerobicPathway]

print(aerobicPathway)

def checkFx(magSet, listOfPathways, threshold):
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
                missing += 1
                if missing == threshold:
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

output = ["\t".join(["MAG_names","Anaerobic", "Aerobic"])]

for mag in misc_list:
    sampleID = "_".join(mag.split("_")[:2])

    presence = ["N","N"]

    for i in range(len(pathwaySet)):
        if checkFx(mag_gene_dict[mag], pathwaySet[i], 5) == True and checkFx(mag_gene_dict[mag], finalPathway, 1):
            presence[i] = "Y"

    output.append("\t".join([mag, "\t".join(presence)]))

with open(f"{basepath}/heatmapGen/heatmapFiles/pathwayCobalamin.txt", 'w') as n:
    n.write("\n".join(output))