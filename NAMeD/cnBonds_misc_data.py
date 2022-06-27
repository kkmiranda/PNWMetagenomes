"""
I used https://iubmb.qmul.ac.uk/enzyme/ to identify the classes of enzymes relevant 
to the generation of ammonia (EC1.4.*, EC3.5.*, EC4.3.1*). Then used those search 
terms in the Orthology database of genome.jp and found all KEGG calls that are genes 
for these classes of enzymes. Now I have 408 enzymes from the KEGG db from these 
classes.
"""

import os

BASEPATH = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/masterdata_17May2021"
FILENAME = "geneFunctionTables/3Jun2022"
kegg_cutoff = 1e-20
category = 'ammonificationHydrolases'

with open(f"{BASEPATH}/heatmapGen/cog_results/{category}.txt","r") as P:
    P = P.readlines()
    gene_dict = {}

    for p in P:
        p = p.split("\t")
        try: 
            gene_dict[p[0]] = " ".join([p[i].strip() for i in range(1,len(p))])
        except:
            print(p)

# Load MAG names
misc_list = {}
with open(f"{BASEPATH}/{FILENAME}/tatoosh_MAG_names.txt",'r') as m:
    m = m.readlines()
    for mag in m:
        misc_list[mag.strip()] = []

mag_gene_dict = {}
with open(f"{BASEPATH}/{FILENAME}/tatoosh_masterdata.txt") as T:
    T = T.readlines()
    for line in T:
        line = line.split()
        if (float(line[-2]) < kegg_cutoff) and ("KOfam" in line[3]): # KEGG e-value
            if line[2] in mag_gene_dict:
                if "!!!" in line[4]:
                    for gene in line[4].split("!!!"):
                        mag_gene_dict[line[2]].append(gene)
                else:
                    mag_gene_dict[line[2]].append(line[4])
            else: # First time mag name appears
                if "!!!" in line[4]:
                    mag_gene_dict[line[2]] = line[4].split("!!!")
                else:
                    mag_gene_dict[line[2]] = [line[4]]

output = []
for mag in misc_list:
    genePresenceCount = 'N'
    for cogGene in gene_dict:
        if cogGene in mag_gene_dict[mag]:
            genePresenceCount = 'Y'
            break
    output.append("\t".join([mag, str(genePresenceCount)]))


header = "MAG_names\tPresence"

output = [header] + output

with open(f"{BASEPATH}/heatmapGen/heatmapFilesPathway/finalAmmonHydrolases.txt", "w") as P:
    P.write("\n".join(output))