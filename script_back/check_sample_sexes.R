library('tidyverse')

meta_df<- read_delim("data/samples.meta.txt", delim = "\t") %>%
  rename(Sample = ebi_sample_acc)

samp_df <- read_delim("data/An1000g_Samples - Sheet3.csv", delim = ",") %>%
  select(Population, Sample)

samp_df %>% left_join(meta_df) %>% 
  write.csv("data/meta_join.csv")
