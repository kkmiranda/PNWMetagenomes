conda activate anvio-7

muscle -in original.fa -out align.fa
trimal -in align.fa -out trim.fa -gt 0.5 -keepheader
anvi-script-reformat-fasta trim.fa --max-percentage-gaps 50 -o g50.fa
iqtree -s g50.fa -nt AUTO -bb 1000 -alrt 1000 -mset LG
anvi-interactive -t *.treefile -p PROFILE.db --manual 
