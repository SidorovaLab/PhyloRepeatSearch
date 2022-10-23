# PhyloRepeatSearch

PhyloRepeatSearch allows the use to run a python script to run a phylogenetic analysis of amino acid repeats across species in comparison to human genes. There are 4 specified hyperparamters that the user can change to generate reports: 

1 - rep_AAs : Give a list of amino acids that appear in the amino acid repeat of interest. For example, ['D','E'], will return reults for any protein sequences that have a continous uniteruppted stretch of Ds and Es. 

2 - gene_ref_accepted_list : Give a list of human genes that also contain the amino acid repeat of interest. For example, ['ATAD2','MYT1','NCL'], will return results for protein othologs across species for those genes.

3 - size : This is the minimum size of the repeat across all species. The default is set to 3, so in the case of reps_AAs of ['D','E'] any occurances of DDD, EEE, DDE, DED, EED, EDE across protein orthologs across species will be returned

4 - gap_allowance : gaps that can interrupt repeated amino acids, default is 0.

For each human protein name given in the gene_ref_accepted_list a corrosponding homology tree is retrieved using the Ensembl rest API (https://rest.ensembl.org/documentation/info/homology_symbol)

Using a recursive algorthim the releative branch length and sequence of each ortholog is retrieved and tabulated. For each otrholog protein sequence retrieved the repeat of length in that sequence is also calculated based on the given hyperparametes.

The protein sequences are then grouped by taxonomy level, giving columns for mean branch length of protein sequences and mean repeat length of protein sequences. This information is stored and exported as phylo_map.csv

Protein sequences for 12 select sepecies (given below) are also stored indiidually by mean branch length of protein sequences and mean repeat length of protein sequences. This information is stored and exported as species_map.csv

saccharomyces_cerevisiae
caenorhabditis_elegans
danio_rerio
xenopus_tropicalis
gallus_gallus
macaca_mulatta
rattus_norvegicus
mus_musculus
bos_taurus
canis_lupus_familiaris
pan_troglodytes
homo_sapiens
