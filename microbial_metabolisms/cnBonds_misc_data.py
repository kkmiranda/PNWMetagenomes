import os

basepath = '/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/masterdata_17May2021'
kegg_cutoff = 1e-20
category = 'ammonificationHydrolases'


with open(f"{basepath}/heatmapGen/cog_results/{category}.txt","r") as P:
    P = P.readlines()
    gene_dict = {}

    for p in P:
        p = p.split("\t")
        try: 
            gene_dict[p[0]] = " ".join([p[i].strip() for i in range(1,len(p))])
        except:
            print(p)

with open(f"{basepath}/geneFunctionTables/goodBins/tatoosh_masterdata.txt") as T:
    T = T.readlines()
    sampleGeneDict = {}

    for t in T:
        t = t.split("\t")
        if "KOfam" in t[3] and float(t[6])<kegg_cutoff:
            gSet = t[4].split("!!!")
            for g in gSet:
                try: 
                    sampleGeneDict[t[2]].add(g)
                except:
                    sampleGeneDict[t[2]] = {g}

output = []
for bin in sampleGeneDict:
    genePresenceCount = 0
    for cogGene in gene_dict:
        sampleID = "_".join(bin.split("_")[:2])
        if cogGene in sampleGeneDict[bin]:
            genePresenceCount+=1
    output.append("\t".join([bin, sampleID, str(genePresenceCount)]))


header = "MAG_names\tenv\tPresence"

output = [header] + output

with open(f"{basepath}/heatmapGen/heatmapFiles/finalAmmonHydrolases.txt", "w") as P:
    P.write("\n".join(output))