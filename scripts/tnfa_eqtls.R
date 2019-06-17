library(TwoSampleMR)
library(tidyverse)

## eqtlgen
tnfceqtl <- 'data/Vosa2018sigciseqtl.chrall.CPRA_b37.tsv.gz' %>%
  read_table2(., comment = '##', col_types = 'cddccdcddddddddccccccc') %>% 
  filter(TRAIT %in% c('TNF', 'LTA', 'TNFRSF1A', 'TNFRSF1B')) %>% 
  mutate(QTL = 'cis_eQTL') %>% 
  write_tsv('input/tnf_ciseqtl.txt')

tnfteqtl <- 'data/Vosa2018sigtranseqtl.chrall.CPRA_b37.tsv.gz' %>%
  read_table2(., comment = '##', col_types = 'cddccdcddddddddccccccc') %>% 
  filter(TRAIT %in% c('TNF', 'LTA', 'TNFRSF1A', 'TNFRSF1B')) %>% 
  mutate(QTL = 'trans_eQTL') %>% 
  write_tsv('input/tnf_transeqtl.txt')

## pqtl
tnfpqtl <- '/Users/sheaandrews/LOAD_minerva/dummy/shea/data/sumstats_munger/output/AholaOlli2018tnfa.chrall.CPRA_b37.tsv.gz' %>%
  read_table2(., comment = '##', col_types = 'cddccdcddddddddccccccc') %>%
  mutate(QTL = 'pqtl') #%>%
  #select(ID, CHROM, POS, REF, ALT, AF, TRAIT, BETA, SE, Z, P, N, DBSNP_ID, OLD_ID) # %>% 
  #write_tsv('/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/0_other/2_DerivedData/tnf_transeqtl.txt')

tnfqtl <- bind_rows(tnfceqtl, tnfteqtl) %>%
  bind_rows(filter(tnfpqtl, P < 1e-4)) %>%
  select(ID, CHROM, POS, REF, ALT, AF, TRAIT, QTL, BETA, SE, Z, P, N, DBSNP_ID, OLD_ID) %>% 
  write_tsv('input/tnf_qtls.txt')




























