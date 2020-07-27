# Core
# Language: R
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 1.1, R 4.0.0
# Dependencies: phyloseq 1.32.0, ape 5.4, psadd 0.1.2, ggplot2 3.3.1, microbiome 1.10.0

PluMA plugin that outputs various statistics on core communities of taxa in a set of samples.

The following are specified in the input TXT file, as tab-delimited keyword-value pairs.

otufile: OTU abundances (CSV)
mapping: Mapping table (CSV)
tree: Phylogenetic tree (CSV)
column: (STRING)
detection: Minimum relative abundance to be counted at all (REAL)
prevalence: Minimum relative abundance to be counted as core

The output prefix will be used for three files:
prefix.members.csv: Core taxa
prefix.abundindex.csv: Percentage of each sample that is a core community
prefix.coverage.csv: Amount of core taxa sequences in each sample

Format of the otufile CSV:
- Rows are OTUs, columns are samples
- Entry (i,j) is the abundance of OTU in sample j.  Can be absolute or relative.

Format of the mapping CSV:
- Rows are samples, columns are metadata
- One single header row, i.e: sample,subject,visit

Format of the tree CSV;
- Rows are OTUs, columns are classifications
- One single header row, i.e. OTU,Kingdom,Phylum,Class,Order,Family,Genus,Species
