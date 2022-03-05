import os

filepath = '/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/masterdata_17May2021/'

Carbohydrates_general = {"COG1653","COG1879","COG3839","COG1175","COG0395","COG1129","COG2211","COG0738","COG1134","COG1682","COG2271","COG5037","COG3822","COG1762","COG1925","COG2213","COG2893","COG3730","COG3732"}
Carbohydrates_pentoses = {"COG1172","COG4214","COG2182","COG3833","COG1455"}
Carboxylic_acids = {"COG4663","COG1593","COG1638","COG4664","COG2358","COG0471","COG4666","COG1301","COG4665","COG0651","COG1620","COG3090","COG1823","COG5037"}
Compatible_solutes = {"COG2113","COG1292","COG4176","COG0591","COG4175","COG1125","COG1174"}

carbo = [
    "COG1653","COG1879","COG3839","COG1175","COG0395","COG1129","COG2211","COG0738","COG1134","COG1682","COG2271","COG5037","COG3822","COG1762","COG1925","COG2213","COG2893","COG3730","COG3732",
    "COG1172","COG4214","COG2182","COG3833","COG1455",
    "COG4663","COG1593","COG1638","COG4664","COG2358","COG0471","COG4666","COG1301","COG4665","COG0651","COG1620","COG3090","COG1823","COG5037",
    "COG2113","COG1292","COG4176","COG0591","COG4175","COG1125","COG1174"
    ]

# Load MAG names
misc_list = {}
with open(f"{filepath}/geneFunctionTables/goodBins/tatoosh_MAG_names.txt",'r') as m:
    m = m.readlines()
    for mag in m:
        misc_list[mag.strip()] = []

# Set up masterdata sheet
mag_gene_dict = {}
with open(f"{filepath}/geneFunctionTables/goodBins/COG_function.txt", 'r') as g:
    g = g.readlines()
    for line in g:
        line = line.split()
        if float(line[-2]) < 1e-50: # COG hits under 1e-50
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

output=[]

for mag in misc_list:
    carbGenSum,carbPentoseSum,carbAcidSum,solute = 0,0,0,0
    for cog in carbo:
        try:
            if cog in mag_gene_dict[mag]:
                if cog in Carbohydrates_general:
                    carbGenSum += 1
                elif cog in Carbohydrates_pentoses:
                    carbPentoseSum += 1
                elif cog in Carboxylic_acids:
                    carbAcidSum += 1
                elif cog in Compatible_solutes:
                    carbAcidSum += 1
        except:
            pass
    
    output.append(f"{mag}\t{carbGenSum}\t{carbPentoseSum}\t{carbAcidSum}\t{carbAcidSum}")
    #print(checking)

header = "MAG_names\tCarbohydrates_general\tCarbohydrates_pentoses\tCarboxylic_acids\tCompatible_solutes"

output.insert(0,header)

with open(f"{filepath}/heatmapGen/heatmapFiles/finalCarbon.txt", "w") as T:
    T.write("\n".join(output))