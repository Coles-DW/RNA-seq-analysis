setwd("D:/OneDrive - Western Sydney University/Donovin PhD/Thesis/Experiments/Chapter 3/Network analysis/SP_host_network/Strong coregulated target plant genes")

#input data
sa <- read.csv ("salicylic acid pathway.csv")
ja <- read.csv ("jasmonic acid pathway.csv")
ep <- read.csv ("ethylene pathway.csv")
pp <- read.csv ("Phenylpropanoid pathway heatmap.csv")

bp <- read.csv ("DEgenes12heatmap.csv")
bns <- read.csv ("degenes24heatmap.csv")
np <- read.csv ("degenes72heatmap.csv")

rd <- read.csv ("resistant drivers heatmap.csv")
sd <- read.csv ("susceptible drivers heatmap.csv")
exp <- read.csv ("Differentialexpression_MR_contRversusContS_v2.csv")
ann <- read.csv ("Egrandis_297_v2.0.annotation_info.csv")

install.packages("dplyr")
library("dplyr")


#merge
sa_bp <- left_join(sa, bp, by="gene_id")
sa_bp_bns <- left_join(sa_bp, bns, by="gene_id")
sa_bp_bns_np <- left_join(sa_bp_bns, np, by="gene_id")

ja_bp <- left_join(ja, bp, by="gene_id")
ja_bp_bns <- left_join(ja_bp, bns, by="gene_id")
ja_bp_bns_np <- left_join(ja_bp_bns, np, by="gene_id")

ep_bp <- left_join(ep, bp, by="gene_id")
ep_bp_bns <- left_join(ep_bp, bns, by="gene_id")
ep_bp_bns_np <- left_join(ep_bp_bns, np, by="gene_id")

pp_bp <- left_join(pp, bp, by="gene_id")
pp_bp_bns <- left_join(pp_bp, bns, by="gene_id")
pp_bp_bns_np <- left_join(pp_bp_bns, np, by="gene_id")


rd_ann <- left_join(rd, ann, by="gene_id")
sd_ann <- left_join(sd, ann, by="gene_id")
exp_ann <- left_join(exp, ann, by="gene_id")

#write files
write.csv(sa_bp_bns_np, "salicylic acid pathway_varieties.csv")
write.csv(ja_bp_bns_np, "jasmonic acid pathway_varieties.csv")
write.csv(ep_bp_bns_np, "ethylene pathway_varieties.csv")
write.csv(pp_bp_bns_np, "phenylpropanoid pathway_varieties.csv")

write.csv(rd_ann, "resistant drivers heat annot.csv")
write.csv(sd_ann, "susceptible drivers heat annot.csv")
write.csv(exp_ann, "annotation of sig genes.csv")
