fq.files <- list.files("/fs/ess/PAS1695/projects/mitox/data/fastqs", full.names = T)

fq.df <- as.data.frame(cbind(fq.files, sampnames = NA)) %>%
  mutate(sampnames = gsub(".*.(FF.*).fastq.gz", "\\1", fq.files))

nfq <- nrow(fq.df)

for(f in 1:nfq){
  fileOut<- paste0("/fs/ess/PAS1695/projects/mitox/scripts/batch/mpa_", fq.df$sampnames[f], ".pbs")
  
  writeLines(c(paste0("#PBS -N mpa_mitox_", fq.df$sampnames[f]),
               "#PBS -A PAS1695",
               "#PBS -l walltime=10:00:00",
               "#PBS -l nodes=1:ppn=28",
               "#PBS -j oe",
               "",
               "cd /fs/ess/PAS1695/projects/mitox/data/metaphlan/",
               "module load python/3.7-2019.10",
               "source activate metaphlan3",
               paste("metaphlan --input_type fastq --read_min_len 50 --add_viruses --nproc 28 --bowtie2db ~/Documents/db/mpa_bowtie2",
                      fq.df$fq.files[f],
                      paste0(fq.df$sampnames[f], ".txt")),
               ""),
             fileOut)
  # close(fileOut)
}

mettab.ls <- lapply(fq.df$sampnames, function(x) read.table(paste0("/fs/ess/PAS1695/projects/mitox/data/metaphlan/", x, ".txt"),
                                                            header = F, sep = "\t", stringsAsFactors = F))
names(mettab.ls) <- fq.df$sampnames
metab.df <- bind_rows(lapply(fq.df$sampnames, function(x) mettab.ls[[x]] %>% mutate(sample = x)))

metab.df.form <- metab.df %>%
  rename("V1" = "Taxonomy",
         "V2" = "TaxNum",
         "V3" = "RelAbun",
         "V4" = "Alternative.Tax")
write.csv(metab.df.form, "/fs/ess/PAS1695/projects/mitox/data/mpa_aggregated_long.csv", row.names = F)
