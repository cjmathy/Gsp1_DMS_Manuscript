# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

df_sets <- read_csv('data/Ras_data/Ras_sca_gof_table.csv', show_col_types=F)
df_annot <- read_csv('scripts/Fig4/toxic_annotations.txt', show_col_types=F)
df_aln <- 
    read_delim('data/Ras_data/HRas_Gsp1_index.txt', delim='\t', show_col_types=F) %>% 
    rename(
        'pos'=Gsp1_res_num,
        'aa_HRas'=HRas,
        'aa_Gsp1'=Gsp1,
        'pos_HRas'=HRas_res_num
    )

df <- 
    left_join(df_sets, df_annot, by=c('pos'='position')) %>%
    left_join(df_aln, by='pos') %>% 
    select(-`C-terminal switch`) %>%  
    rename('PTM sites in Gsp1'=PTMs) %>% 
    mutate(
        `GTPase regions`=ifelse(`GTPase regions`==0, F, T),
        `Contacts Nucleotide`=ifelse(`Contacts Nucleotide`==0, F, T),
        `Distal sites affecting switching`=ifelse(`Distal sites affecting switching`==0, F, T),
        `PTM sites in Gsp1`=ifelse(`PTM sites in Gsp1`==0, F, T),
        `Regulator Interface`=ifelse(`Regulator Interface`==0, F, T),
    ) %>% 
    select(-gof_all) %>% 
    rename('gof'=gof_2plus) %>% 
    pivot_longer(cols=-c(pos,tox,gof,sca,aln_num,aa_HRas,aa_Gsp1,pos_HRas), names_to='category', values_to='in_region') %>% 
    mutate(category=factor(category, levels=c(
        'GTPase regions',
        'Contacts Nucleotide',
        'Distal sites affecting switching',
        'PTM sites in Gsp1',
        'Regulator Interface'
    ))) %>% 
    mutate(res_Gsp1 = paste0(aa_Gsp1,pos), res_HRas = paste0(aa_HRas,pos_HRas)) %>% 
    select(res_Gsp1, res_HRas, tox, gof, sca, category, in_region) %>% 
    nest(data = c(res_Gsp1, res_HRas, tox, gof, sca)) %>% 
    filter(in_region) %>% 
    select(-in_region) %>% 
    mutate(
        pasted = map(data, function(d) {
                tibble(
                    'tox'=filter(d, tox) %>% pull(res_Gsp1) %>% paste(collapse = ', '),
                    'gof'=filter(d, gof) %>% pull(res_HRas) %>% paste(collapse = ', '),
                    'sca'=filter(d, sca) %>% pull(res_HRas) %>% paste(collapse = ', ')
                )
            }
        )
    ) %>% 
    unnest(pasted) %>% 
    select(-data) %>% 
    arrange(category) %>% 
    rename(
        'Gsp1 residue category (Figure 4)' = category,
        'Toxic/GOF positions in Gsp1 (10 or more toxic mutations)' = tox,
        'Activating positions in HRas (2 or more activating mutations)' = gof,
        'Sector positions in HRas (from SCA analysis)' = sca,
    )

df

write_csv(df, 'figures/Tables/TableS1_Gsp1_HRas_annotations.csv')




## Show the alignment for HRas - Gsp1

