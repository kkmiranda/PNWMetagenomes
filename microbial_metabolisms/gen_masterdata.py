import os 
import re

"""
I take in the gene function file, add in:
1. split name
2. bin name
3. tax_id
"""

## Would have to update the basepaths to your configuration
gene_basepath = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/gene_splits"
bin_basepath = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/bins_n_splits"
func_basepath = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/gene_function"
taxon_basepath = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/tax_ID"
output_basepath = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/masterdata_17May2021"

def gen_masterdata(Sample):
    # def gen_masterdata(Sample):
    gcall_file = f"{Sample}_gene_splits.txt" #
    bin_file = f"{Sample}_concoct.txt" # 
    func_file = f"{Sample}_functions.txt"
    tax_file = f"{Sample}_taxID.txt"
    output_file = f"rawMasterData/{Sample}_masterdata.txt"

    #######
    """
    Get gene calls: Gene call file
    From: anvi-export-table CONTIGS.db --table genes_in_splits

    Gene_call   split  gene_length
    """
    with open(f"{gene_basepath}/{gcall_file}","r") as g:
        g = g.readlines()
        g.pop(0)

        gene_split = {} 
        for gene in g:
            gene = gene.split("\t")
            split_name = gene[0]
            gene_id = gene[1]
            # start = int(gene[2])
            # stop = int(gene[3])
            # gene_length = abs(stop - start)
            gene_split[gene_id] = split_name # keyed by split name
            
    #OUTPUT DICT: gene_split
    #######

    #######
    """
    Get TAXONOMY
    From: anvi-estimate-scg-taxonomy

    Bin_name    taxonomic_ID
    """
    with open(f"{taxon_basepath}/{tax_file}","r") as t:
        t = t.readlines()
        t.pop(0)

        bins_taxID = {}
        for bin in t:
            bin = bin.split("\t")
            bin_name = bin[0]
            # tot_scg = bin[1]
            # sup_scg = bin[2]
            
            taxID = [x.strip() for x in bin[3:] if x!='None'] #clean up the taxon names
            if len(taxID) == 0: #in the case that the bin just isn't resolving itself
                taxID = "NA"

            taxID = "_".join(taxID)
            taxID = taxID.replace(" ","_")
            # bin_tuple = (
            #     tot_scg,
            #     sup_scg,
            #     taxID
            # )
            # bins_taxID[bin_name] = bin_tuple
            bins_taxID[bin_name] = taxID

    #######

    #######
    """
    Get corresponding split and bin: Bins_n_split file
    From: anvi-export-collection

    Split_name      Bin_name
    """
    with open(f"{bin_basepath}/{bin_file}","r") as b:
        b = b.readlines()

        split_bin = {}
        for split in b:
            split = split.split("\t")
            split_name = split[0]
            bin_name = split[1].strip()
            split_bin[split_name] = bin_name

    # print("Split_Bins", len(split_bin))
    #OUTPUT DICT: split_bins
    #######

    #######

    with open(f"{func_basepath}/{func_file}", 'r') as f:
        f = f.readlines()
        f.pop(0)

        f = [(line.strip()).split("\t") for line in f]

    output_data = []
    for gene in f:
        try:
            gene_id = gene.pop(0)
            split_name = gene_split[gene_id]
            bin_name = split_bin[split_name]
            try:
                tax_id = bins_taxID[bin_name]
            except:
                tax_id = "Unknown_taxon"
            end_list = [gene_id, split_name, bin_name, "\t".join(gene), tax_id]

            output_data.append("\t".join(end_list))
        except:
            pass

    with open(f"{output_basepath}/{output_file}",'w') as t:
        t.write("\n".join(output_data))

if __name__=="__main__":
    for sample in ["MG1","MG2","MG3","MG4","MG5","MG6","MG7","MG8"]:
        gen_masterdata(sample)
