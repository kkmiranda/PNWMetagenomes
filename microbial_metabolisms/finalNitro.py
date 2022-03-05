import os
basepath = os.getcwd()

filepath = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/masterdata_17May2021"

Dissimilatory_nitrate_reduction = {"K00370","K00371","K00374","K02567","K02568","K00362","K00363","K03385","K15876"}
Assimilatory_nitrate_reduction = {"K00367","K10534","K00372","K00360","K17877","K00366"}
Denitrification = {"K00370","K00371","K00374","K02567","K02568","K00368","K15864","K04561","K02305","K00376"}
Nitrogen_fixation = {"K02586","K02591","K02588","K00531","K22896","K22897","K22898","K22899"}
Nitrification = {"K10944","K10945","K10946","K10535","K00370","K00371"}
Annamox = {"K10944", "K00368","K15864","K20935","K20932","K20933","K20934"}
Urease = {"K01427","K01428","K01429","K01430","K14048"}

nitro = ["K00370","K00371","K00374","K02567","K02568","K00362","K00363","K03385","K15876","K00367","K10534","K00372","K00360","K17877","K00366","K00370","K00371","K00374","K02567","K02568","K00368","K15864","K04561","K02305","K00376","K02586","K02591","K02588","K00531","K22896","K22897","K22898","K22899","K10944","K10945","K10946","K10535","K00370","K00371","K10944", "K00368","K15864","K20935","K20932","K20933","K20934","K01427","K01428","K01429","K01430","K14048"]


# Load MAG names
misc_list = {}
with open(f"{filepath}/geneFunctionTables/goodBins/tatoosh_MAG_names.txt",'r') as m:
    m = m.readlines()
    for mag in m:
        misc_list[mag.strip()] = []

# Set up masterdata sheet
mag_gene_dict = {}
with open(f"{filepath}/geneFunctionTables/goodBins/tatoosh_masterdata.txt", 'r') as g:
    g = g.readlines()
    for line in g:
        line = line.split()
        if float(line[-2]) < 1e-20: # KEGG hits under 1e-20
            if line[2] in mag_gene_dict:
                if "!!!" in line[4]:
                    for kegg in line[4].split("!!!"):
                        mag_gene_dict[line[2]].append(kegg)
                else:
                    mag_gene_dict[line[2]].append(line[4])
            else: # First time mag name appears
                if "!!!" in line[4]:
                    mag_gene_dict[line[2]] = line[4].split("!!!")
                else:
                    mag_gene_dict[line[2]] = [line[4]]

for kegg in nitro:
    for mag in misc_list:
        try: 
            if kegg in mag_gene_dict[mag]:
                misc_list[mag].append("Y")
            else:
                misc_list[mag].append("N")
        except:
            misc_list[mag].append("N")

output = [key + "\t" + "\t".join(misc_list[key]) for key in misc_list]

header = [
    "MAG_names",
    "Dissimilatory nitrate reduction_K00370","Dissimilatory nitrate reduction_K00371","Dissimilatory nitrate reduction_K00374","Dissimilatory nitrate reduction_K02567","Dissimilatory nitrate reduction_K02568","Dissimilatory nitrate reduction_K00362","Dissimilatory nitrate reduction_K00363","Dissimilatory nitrate reduction_K03385","Dissimilatory nitrate reduction_K15876",
    "Assimilatory nitrate reduction_K00367","Assimilatory nitrate reduction_K10534","Assimilatory nitrate reduction_K00372","Assimilatory nitrate reduction_K00360","Assimilatory nitrate reduction_K17877","Assimilatory nitrate reduction_K00366",
    "Denitrification_K00370","Denitrification_K00371","Denitrification_K00374","Denitrification_K02567","Denitrification_K02568","Denitrification_K00368","Denitrification_K15864","Denitrification_K04561","Denitrification_K02305","Denitrification_K00376",
    "Nitrogen fixation_K02586","Nitrogen fixation_K02591","Nitrogen fixation_K02588","Nitrogen fixation_K00531","Nitrogen fixation_K22896","Nitrogen fixation_K22897","Nitrogen fixation_K22898","Nitrogen fixation_K22899",
    "Nitrification_K10944","Nitrification_K10945","Nitrification_K10946","Nitrification_K10535","Nitrification_K00370","Nitrification_K00371",
    "Annamox_K10944","Annamox_K00368","Annamox_K15864","Annamox_K20935","Annamox_K20932","Annamox_K20933","Annamox_K20934",
    "Urease_K01427","Urease_K01428","Urease_K01429","Urease_K01430","Urease_K14048"
    ]

output.insert(0,"\t".join(header))
# Currently have to manually rename headers because certain pathways have the same genes

with open(f"{filepath}/heatmapGen/heatmapFiles/finalNitrogen.txt", "w") as T:
    T.write("\n".join(output))
           
