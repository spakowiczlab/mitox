library(tidyverse)

fq.files <- list.files("/fs/ess/PAS1695/projects/mitox/data/2022-08-08_fastqs", full.names = T)

fq.df <- as.data.frame(cbind(fq.files, sampnames = NA)) %>%
  mutate(sampnames = gsub(".*.(FF.*).fastq.gz", "\\1", fq.files))

nfq <- nrow(fq.df)

for(f in 1:nfq){
  fileOut<- paste0("/fs/ess/PAS1695/projects/mitox/scripts/batch3/humann3_", fq.df$sampnames[f], ".pbs")

  writeLines(c(paste0("#PBS -N humann3_mitox_", fq.df$sampnames[f]),
               "#PBS -A PAS1695",
               "#PBS -l walltime=10:00:00",
               "#PBS -l nodes=1:ppn=10",
               "#PBS -j oe",
               "",
               "cd /fs/ess/PAS1695/projects/mitox/data/2022-08-08_humann3/",
               "module load python/3.7-2019.10",
               "source activate humann3",
               paste0("humann -i ", fq.df$fq.files[f]," -o ", fq.df$sampnames[f], " --threads 10 --input-format fastq.gz",
                      " --nucleotide-database /fs/ess/PAS1695/db/chocophlan/chocophlan --protein-database /fs/ess/PAS1695/db/chocophlan/uniref"),
               ""),
             fileOut)
  # close(fileOut)
}

mettab.ls <- lapply(fq.df$sampnames, function(x) read.table(paste0("/fs/ess/PAS1695/projects/mitox/data/2022-08-08_humann3/",
                                                                   x, "/",x, "_humann_temp/", x, "_metaphlan_bugs_list.tsv"),
                                                            header = F, sep = "\t", stringsAsFactors = F) %>%
                      mutate(V2 = as.character(V2)))
names(mettab.ls) <- fq.df$sampnames
metab.df <- bind_rows(lapply(fq.df$sampnames, function(x) mettab.ls[[x]] %>% mutate(sample = x)))



metab.df.form <- metab.df %>%
  rename("Taxonomy" = "V1",
         "V2" = "TaxNum" = "V2",
         "RelAbun" = "V3",
         "Alternative.Tax" = "V4")
write.csv(metab.df.form, "/fs/ess/PAS1695/projects/mitox/data/2022-08-08_mpa_aggregated_long.csv", row.names = F)


humann.ls <- lapply(fq.df$sampnames, function(x) read.table(paste0("/fs/ess/PAS1695/projects/mitox/data/2022-08-08_humann3/",
                                                                   x, "/",x, "_pathabundance.tsv"),
                                                            header = F, sep = "\t", stringsAsFactors = F) %>%
                      mutate(sample = x))

humann3.df <- humann.ls %>%
  bind_rows() %>%
  rename("Pathway" = "V1",
         "PathAbun" = "V2")
write.csv(humann3.df, "/fs/ess/PAS1695/projects/mitox/data/2022-08-08_humann3-aggreagated.csv", row.names = F)
