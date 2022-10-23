import requests
import numpy as np
import pandas as pd
import time
import warnings
 
rep_AAs = ['D','E'] # Put amino acids present in amino acid repeat of intrest in list format
gene_ref_accepted_list = ['ATAD2','MYT1','NCL'] # Select subset of human genes which contain the amino acid repeat of intrest
size = 3 # minimum repeat size across all species, should be more than or equal to 3
gap_allowance = 0 # gaps inbetween repeated amino acids, default is 0

t = time.time()
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

def time_run(t):
    elapsed = time.time() - t
    elapsed_min = elapsed / 60
    elapsed_min = "{:.2f}".format(elapsed_min)
    return (str(elapsed_min) + ' min')

phylo_map = pd.DataFrame(
    columns = (pd.MultiIndex.from_product([
        gene_ref_accepted_list,
        ['branch_length','repeat_length']])))

species_map = pd.DataFrame(
    columns = (pd.MultiIndex.from_product([
        gene_ref_accepted_list,
        ['repeat_length']])))

gene_num = -1

for gene in gene_ref_accepted_list:
    
    print('\nOrthologs of ' + gene + ': ' + time_run(t))
    
    server = "https://rest.ensembl.org"
    ext = "/homology/symbol/human/" + gene + "?format=condensed"
     
    r = requests.get(server+ext, 
                     headers={ "Content-Type" : "application/json"})
     
    if not r.ok:
        print('server error')
        
    else:  
        
        gene_num += 1
        decoded = r.json()
        homologs = decoded["data"][0]
        
        hdf = pd.DataFrame(homologs["homologies"])
        
        hdf.drop(hdf.index[hdf['type'] == 'within_species_paralog'],
                 inplace = True)
        
        hdf.drop(hdf.index[hdf['type'] == 'other_paralog'],
                 inplace = True)    
        
        ext = ("/genetree/member/symbol/homo_sapiens/" 
            + gene 
            + "?sequence=none;content-type=application/json")
         
        r = requests.get(server+ext, 
                         headers={ "Content-Type" : "application/json"})
         
        if not r.ok:
          print('server error')
          hdf = hdf.drop(hdf.index[gene_num])
         
        decoded = r.json()
        
        tree = decoded["tree"]
        tax = tree['taxonomy']
        
        global treeDF
        treeDF = pd.DataFrame(
                {'species': tax['scientific_name'],
                 'branch_length': tree['branch_length']},
                 index = [tax['id']])
        
        print('Calculating relation tree: ' + time_run(t))
        
        def append_df(tree):
            global treeDF
            tree = tree['children']
            
            for child in tree:
                
                tax = child['taxonomy']
                name = tax['scientific_name']
                name = name.replace(" ", "_")
                name = name.replace("_strain_S288C", "")
                name = name.replace("_strain_N2", "")
                name = name.replace("_reference_(cl57bl6)_strain", "")
    
                stip_name = name.lower()
                
                treeDF = treeDF.append(
                    pd.DataFrame(
                        {'species': stip_name,
                         'branch_length': child['branch_length']},
                        index = [tax['id']]))
                
                if 'children' in child:
                    append_df(child)
                
        append_df(tree)
        
        treeDF = treeDF.sort_values('branch_length', 
            ascending=False).drop_duplicates(['species'])
        
        treeDF = treeDF.reset_index()
        
        hdf = pd.merge(hdf, treeDF)
        hdf = hdf.drop('index', 1)
        hdf = hdf.fillna(0)
        
        spec_ids = hdf[['species']].values.tolist()
        spec_ids = [i[0] for i in spec_ids]
        
        prot_ids = hdf[['protein_id']].values.tolist()
        prot_ids = [i[0] for i in prot_ids]
        
        print('\nRetrieving sequences: ')
        row = -1
            
        for id in prot_ids:
            
            print(spec_ids[row] + " (" + gene + ")")
            row += 1
            ext = "/sequence/id/" + id + "?type=protein"
            r = requests.get(server+ext, 
                             headers={ "Content-Type" : "application/json"})
             
            if not r.ok:
                print('server error')
                hdf = hdf.drop(hdf.index[row])
             
            else:
                decoded = r.json()
                seq = (repr(decoded["seq"]))
                hdf.loc[row, 'seq'] = seq 
        
        seq_list = hdf[['seq']].values.tolist()
        seq_list = [i[0] for i in seq_list]
        
        print('\nRepeat search: ' + time_run(t))
        
        row = -1
        
        for seq in seq_list:
            
            row += 1
            domain_len = 0 ; res_num_cont = 0 ; gap_num = 0
            rep_len_cont = 0; prev_count = 0; rep_len_cont_final = 0
    
            for AA_pos in range(len(seq)):
                domain_len += 1 ; rep_len_cont += 1
                
                if seq[(domain_len-1) : domain_len] in rep_AAs :
                    res_num_cont += 1 ; gap_num = 0
    
                else:
                    gap_num += 1
    
                    if gap_num > gap_allowance:
    
                        if res_num_cont >= size:
                            
                            if rep_len_cont > prev_count:
                            
                                prev_count = rep_len_cont
                                rep_len_cont_final = rep_len_cont - 1
    
                        res_num_cont = 0; rep_len_cont = 0
            
            print(spec_ids[row] + " (" + gene + "): " + str (rep_len_cont_final))
            hdf.loc[row, 'repeat_length'] = rep_len_cont_final
                
        hdf = hdf.sort_values('repeat_length', 
                              ascending=False).drop_duplicates(['species'])
        
        hdf = hdf.replace(-1,0)
        hdf = hdf.replace(1,0)
            
        hdf_stats_phylo = hdf[['taxonomy_level','branch_length','repeat_length']]
        hdf_stats_phylo = hdf_stats_phylo.groupby('taxonomy_level').mean()
        hdf_stats_species = hdf[['species','repeat_length']]
        hdf_stats_species = hdf_stats_species[
                (hdf_stats_species["species"] == 'saccharomyces_cerevisiae') | 
                (hdf_stats_species["species"] == 'caenorhabditis_elegans') | 
                (hdf_stats_species["species"] == 'danio_rerio') | 
                (hdf_stats_species["species"] == 'xenopus_tropicalis') | 
                (hdf_stats_species["species"] == 'gallus_gallus') | 
                (hdf_stats_species["species"] == 'macaca_mulatta') | 
                (hdf_stats_species["species"] == 'rattus_norvegicus') | 
                (hdf_stats_species["species"] == 'mus_musculus') | 
                (hdf_stats_species["species"] == 'bos_taurus') | 
                (hdf_stats_species["species"] == 'canis_lupus_familiaris') | 
                (hdf_stats_species["species"] == 'pan_troglodytes') | 
                (hdf_stats_species["species"] == 'homo_sapiens')
                ]
        
        hdf_stats_phylo.columns = pd.MultiIndex.from_product([[gene], hdf_stats_phylo.columns])
        phylo_map = phylo_map.append(hdf_stats_phylo)
        phylo_map = phylo_map.groupby(level=0).sum()
        
        hdf_stats_species = hdf_stats_species.set_index('species')
        hdf_stats_species.columns = pd.MultiIndex.from_product([[gene], hdf_stats_species.columns])
        species_map = species_map.append(hdf_stats_species)
        species_map = species_map.groupby(level=0).sum()
        
        for index, row in phylo_map.iterrows():
            
            for i in range(int(len(row)/2)):
                
                if row[2*i] == 0 and row[(2*i)+1] == 0:
                    row[2*i] = np.nan
                    row[(2*i)+1] = np.nan
                    
species_map.to_csv('species_map.csv')
phylo_map.to_csv('phylo_map.csv')  