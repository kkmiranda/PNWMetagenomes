# The diversity and functional capacity of microbes associated with coastal macrophytes
### <b>Authors</b>: Khashiff K. Miranda (corresponding), Brooke L. Weigel, Emily C. Fogarty, Iva A. Veseli, Anne E. Giblin, A. Murat Eren, Catherine A. Pfister

<br>
This folder contains files and workflows used in describing novel microbial communities living on the surfaces and in the rhizomes of marine coastal phototrophs in the Pacific NorthEast Ocean. Sampling was conducted on Tatoosh Island (Makah Tribe, WA, USA) and at West Falmouth Harbor (MA, USA). 

- Start with <a href='./01_Raw_Sequences_to_MAGs.md'>From Raw Sequences to MAGs</a>
- <a href='./fig1_DissolvedO2/README.md'>Fig. 1: Dissolved Oxygen Graph</a>
- <a href='./fig2_taxonomy/README.Rmd'>Fig. 2: Taxonomy</a>
- <a href='./NAMeD/README.md'>Fig. 3: Microbial Metabolism Heatmap using NAMeD</a>
- <a href='./fig4nifH_phylogeny/README.md'>Fig. 4: Identifying nifH genes and generating a phylogenetic tree</a>


### Database Files available at:
Figshare: 10.6084/m9.figshare.20152949
    - If you want to follow all our steps from scratch (Raw sequences to figures), you'll have to download all the contig and profile database files from this repo. Make new directories named 03_CONTIGS/ and 05_ANVIO_PROFILES/ rooted into this directory and store files in there accordingly. 
    - If you want to just mess around with NAMeD and generate the figures (particularly Fig 3), just download `tatoosh_masterdata.txt` and `tatoosh_MAG_names.txt`. I'd save these two files into the folder named assets/

### Workflow credits
Organisation of this workflow was inspired by M. Eren's countless <a href='https://anvio.org/'>anvi'o tutorials</a> and B. Weigel's <a href='https://github.com/brookeweigel/Kelp_associated_bacterial_genomes'>Kelp Metagenome Workflow</a>

## Contact
For any queries, feel free to contact khashiff.m@gmail.com or cpfister@uchicago.edu