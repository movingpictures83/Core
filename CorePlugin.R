library(microbiome)
library(ggplot2)
#library(phyloseq)
library(ape)
library(psadd)

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1]; 
   # Need to get the three files
   otu.path <<- parameters["otufile", 2]
   tree.path <<- parameters["tree", 2]
   map.path <<- parameters["mapping", 2]
   detection <<- parameters["detection", 2]
   preval <<- parameters["prevalence", 2]
   #HMP <<- import_qiime(otu.path, map.path, tree.path, parseFunction = parse_taxonomy_qiime)
}
run <- function() {
   #samples.to.keep <<- sample_sums(HMP) >= 1000
   #HMP <<- prune_samples(samples.to.keep, HMP)
   #HMP <<- filter_taxa(HMP, function(x) sum(x >3) > (0.01*length(x)), TRUE)
   physeq <<- read_csv2phyloseq(otu.file=otu.path, taxonomy.file=tree.path, metadata.file=map.path)
   mytree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
   physeq <<- merge_phyloseq(physeq, mytree)
}
output <- function(outputfile) {
  #height = 10*300); #,)
  #result <<- PCoA(physeq)
  abundindex <- core_abundance(physeq, detection, preval)
  members <- core_members(physeq, detection, preval)
  cover <- coverage(physeq, threshold = preval) 
  pdf(paste(outputfile,"pdf",sep="."))
  y <- plot_core(physeq, prevalences=seq(detection, 1, 0.1), min.prevalence=prevalence)
  #y <- plot_sparsity(p0)
  #print(str(y))
  #print(str(y$data))
  write.csv(members, paste(outputfile, "members", "csv", sep="."))
  write.csv(abundindex, paste(outputfile, "abundindex", "csv", sep="."))
  write.csv(cover, paste(outputfile, "coverage", "csv", sep="."))
  print(y)#plot_bar(HMP, x="Description", fill=diffcol))
  dev.off()

}
#input("plugins/Bar/example/parameters.txt")
#run()
#output("plugins/Bar/example/yes.pdf")

