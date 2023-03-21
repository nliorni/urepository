#!/usr/bin/env Rscript

suppressMessages(library(dplyr))


args <- commandArgs(T)
if(args[1] == '-h'){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt))
    print("Usage:")
    print("MEQC.R [thresholds] [output dir] [samtools input] [bedtools input]")
    print("Checks complete. Version 0.1")
    stop()
}

if (is.na(args[1])){
  print("No thresholds provided, using default ones: 1,2,3,5,10,20,30,50,100,500,1000,5000")
  thresholds <- c(1,2,3,5,10,20,30,50,100,500,1000,5000)
}else{
  thresholds <- as.numeric(unlist(strsplit(args[1], split=",")))
}
path_out <- args[2]
samstats <- read.csv(args[3], sep='\t', header=FALSE)
positions <- as.numeric(length(samstats$V1))
sink(path_out)
writeLines("Sensitivity\n\nCov\tTargetFraction")
for (c in thresholds){
  writeLines(paste(c, length( which( samstats$V3 > c ) )/positions, sep = '\t'))
}

writeLines("\nStats\n")
writeLines(paste("Mean Cov", mean(samstats$V3), sep = '\t'))
writeLines(paste("StDev",sd(samstats$V3), sep = '\t\t'))
writeLines("Quartiles")
print(quantile(samstats$V3))

bedstats <- read.csv(args[4], sep = '\t', header=FALSE)
bedstats[5] <- NULL

targetstats <- bedstats %>%
  group_by(V1, V2, V3, V4) %>%
  summarise_all(funs(mean, sd))

targetstats <- targetstats %>% mutate_if(is.numeric, funs(round(., 3)))

targetstats$Size <- targetstats$V3 - targetstats$V2
writeLines('\nChr\tStart\tEnd\tFeatureName\tMeanCov\tStDev\tSize\n')
write.table(targetstats, file = path_out, quote= FALSE, append = TRUE, sep= '\t', row.names=FALSE, col.names=FALSE)
sink()




