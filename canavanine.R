## Canavanine Microbiome

# Set working directory
setwd("~/Desktop/2021-0210/Canavanine-R-Tutorial/")

# Install and load some packages
if (!require("gplots")) install.packages("gplots")
if (!require("stringr")) install.packages("stringr")
library("gplots")
library("stringr")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("ALDEx2")) BiocManager::install("ALDEx2")
library("ALDEx2")

if (!require("RColorBrewer")) install.packages("RColorBrewer")
library("RColorBrewer")

# Make color palettes
RColorBrewer::display.brewer.all()
cols <- brewer.pal(9, "Set1")
cols2 <- brewer.pal(8, "Set2")
cols3 <- brewer.pal(12, "Set3")
cols13 = c(cols[c(1:5,7:9)],"black",cols2[c(1,3,5,6)]) 
cols14 = c(cols,"black",cols2[c(1,3,5,6)]) 
gcol = colorRampPalette(c("black", "white"))(14)[3:12]


# Sample metadata
metadata_cnv = read.table("Input_data/metadata_cana.txt", header = F, as.is = T, sep = "\t")
colnames(metadata_cnv) = c("Sample_ID","Chemical","Conc","Group")


# beta diversity PCoA
wuf_pcoa = read.table("Input_data/weighted_unifrac_pcoa_ordination/ordination.txt", skip = 9, nrows = 16, sep = "\t")
wuf_pcoa_explain = read.table("Input_data/weighted_unifrac_pcoa_ordination/ordination.txt", skip = 4, nrows = 1, sep = "\t")
uuf_pcoa = read.table("Input_data/unweighted_unifrac_pcoa_ordination/ordination.txt", skip = 9, nrows = 16, sep = "\t")
uuf_pcoa_explain = read.table("Input_data/unweighted_unifrac_pcoa_ordination/ordination.txt", skip = 4, nrows = 1, sep = "\t")
wuf_pcoa = wuf_pcoa[c(1:4,9:16,5:8),]
uuf_pcoa = uuf_pcoa[c(1:4,9:16,5:8),]

plot(x = wuf_pcoa$V2, y = wuf_pcoa$V3)
plot(x = uuf_pcoa$V2, y = uuf_pcoa$V3)

# the paste function is to select the axcel of ecach PCo precentage eplaned.

file_name = "Output_data/Canavanine_WUF_UUF_PCoA.pdf"
pdf(file_name, width=9.0, height=4.2)
layout(matrix(c(1,1,1,2,2,2,3),nrow = 1))
plot(x = wuf_pcoa$V2, y = wuf_pcoa$V3, xlab = paste("PCo1 (",round(wuf_pcoa_explain$V1,3)*100,"%)",sep = ""),
     ylab = paste("PCo2 (",round(wuf_pcoa_explain$V2,3)*100,"%)",sep = ""), pch = 21,
     col = c(adjustcolor(cols13[rep(5,8)], alpha = 0.9), adjustcolor(cols13[rep(9,8)],alpha.f = 0.9)), 
     bg = c(rep(NA,4), adjustcolor(cols13[rep(5,4)], alpha = 0.2),
            adjustcolor(cols13[rep(5,4)], alpha = 0.4), adjustcolor(cols13[rep(5,4)], alpha = 0.8)),
     main = "Weighted UniFrac PCoA", cex =2)
plot(x = uuf_pcoa$V2, y = uuf_pcoa$V3, xlab = paste("PCo1 (",round(uuf_pcoa_explain$V1,3)*100,"%)",sep = ""),
     ylab = paste("PCo2 (",round(uuf_pcoa_explain$V2,3)*100,"%)",sep = ""), pch = 21,
     col = c(adjustcolor(cols13[rep(5,8)], alpha = 0.9), adjustcolor(cols13[rep(9,8)],alpha.f = 0.9)), 
     bg = c(rep(NA,4), adjustcolor(cols13[rep(5,4)], alpha = 0.2),
            adjustcolor(cols13[rep(5,4)], alpha = 0.4), adjustcolor(cols13[rep(5,4)], alpha = 0.8)),
     main = "Unweighted UniFrac PCoA", cex =2)

#in the next step we add le legend to the graph and made a lagend layer

plot(x=1, y=1, ann = F, xaxt="n", yaxt="n", bty = "n", ty = "n")
legend("right", legend = c("Control","Low","Middle","High"), xpd = T, pch = 21, col = adjustcolor(cols13[c(5,5,9,9)],alpha.f = 0.9), 
       pt.bg = c(NA, adjustcolor(cols13[5], alpha = 0.2), adjustcolor(cols13[5], alpha = 0.4), adjustcolor(cols13[5], alpha = 0.8)))
dev.off()


# Alpha diveristy data are obtained from the qiime2 wivw and safed as .tsv
shannon_meta = read.table("Input_data/shannon_metadata.tsv", header = T, as.is = T,sep = "\t")
obs_otu_meta = read.table("Input_data/obsOTU_metadata.tsv", header = T, as.is = T,sep = "\t")
faith_meta = read.table("Input_data/faithPD_metadata.tsv", header = T, as.is = T,sep = "\t")

boxplot(observed_otus ~ Conc, data = obs_otu_meta)
boxplot(shannon ~ Conc, data = shannon_meta)
boxplot(faith_pd ~ Conc, data = faith_meta)


#str_sub function dfrom the sting pakage is for the costomizing the gopu names in the x lab
  
file_name = "Output_data/Canavanine_Alpha_div.pdf"
pdf(file_name, width=10, height=4.0)
par(mfrow = c(1,3))
boxplot(observed_otus ~ Conc, data = obs_otu_meta,  names = str_sub(sort(unique(obs_otu_meta$Conc)), start = 4),
        col =  adjustcolor(cols13[5], alpha = 0.8), pch = 21, cex = 0.1, 
        main = "Observed ASVs", xlab = "", ylab = "Observed ASVs", las = 1)
par(new = T)
plot(x = rep(c(1,4,2,3),each = 4), y = obs_otu_meta$observed_otus, 
     xaxt = "n",yaxt ="n",bty = "n", xlab = "", ylab = "", pch = 21, col = 1, bg = "white", xlim = c(0.5,4.5), cex = 0.8)
text(x = 3:4, y = c(max(obs_otu_meta$observed_otus[13:16]),max(obs_otu_meta$observed_otus[5:8])) + 50,"*", cex = 1.2)

boxplot(shannon ~ Conc, data = shannon_meta,  names = str_sub(sort(unique(obs_otu_meta$Conc)), start = 4),
        col =  adjustcolor(cols13[5], alpha = 0.8), pch = 21, cex = 0.1, 
        main = "Shannon", xlab = "", ylab = "shannon index", las = 1)
par(new = T)
plot(x = rep(c(1,4,2,3),each = 4), y = shannon_meta$shannon, 
     xaxt = "n",yaxt ="n",bty = "n", xlab = "", ylab = "", pch = 21, col = 1, bg = "white", xlim = c(0.5,4.5), cex = 0.8)

boxplot(faith_pd ~ Conc, data = faith_meta,  names = str_sub(sort(unique(obs_otu_meta$Conc)), start = 4),
        col =  adjustcolor(cols13[5], alpha = 0.8), pch = 21, cex = 0.1, 
        main = "Faith's PD", xlab = "", ylab = "faith PD index", las = 1)
par(new = T)
plot(x = rep(c(1,4,2,3),each = 4), y = faith_meta$faith_pd, 
     xaxt = "n",yaxt ="n",bty = "n", xlab = "", ylab = "", pch = 21, col = 1, bg = "white", xlim = c(0.5,4.5), cex = 0.8)
text(x = 3:4, y = c(max(faith_meta$faith_pd[13:16]),max(faith_meta$faith_pd[5:8])) + 2,"*", cex = 1.2)
dev.off()



# Taxonomy, Family level, at fisrt we clean the name of each taxa by selcting the fist row  [1,] and to remove the ; or ) , 
#apply function is fo removing the 0 from the cilumns 
taxon_L5 = read.csv("Input_data/level-5_s132.csv", header = T, as.is = T)
taxon_L5_head = read.csv("Input_data/level-5_s132.csv", header = F, as.is = T)[1,]
colnames(taxon_L5) = taxon_L5_head
numcol = ncol(taxon_L5)
taxon_L5_cnv = cbind(taxon_L5[,c(1, (numcol-2):numcol)], 
                     taxon_L5[,2:(numcol-3)][,apply(taxon_L5[,2:(numcol-3)],2,max) > 0])

# Transform Read count --> Relative abundance (%)
rate_L5_cnv = cbind(taxon_L5_cnv[,1:4], taxon_L5_cnv[,-(1:4)]/rowSums(taxon_L5_cnv[,-(1:4)])*100)
# Filter out low abundant bacteria (less reliable)
rate_L5_cnv_filt = cbind(rate_L5_cnv[,1:4],rate_L5_cnv[,-(1:4)][,colMeans(rate_L5_cnv[,-(1:4)]) > 0.05])

# Differential abundance test, Welch's t-test corrected by Benjamini-Hochberg method
# we use the for function to file up the whole matrix tha we made for the p value and log2  in the t_result  then we adjust the p value to q vlue  using the p.ajust function.
t_result = matrix(NA, ncol = 3, nrow = ncol(rate_L5_cnv_filt)-4)
for (i in 1:(ncol(rate_L5_cnv_filt)-4)) {
  t_result[i,1] = t.test(rate_L5_cnv_filt[1:4,i+4], rate_L5_cnv_filt[5:8,i+4], var.equal = F, paired = F)$p.value
  t_result[i,3] = log2(mean(rate_L5_cnv_filt[5:8,i+4])/mean(rate_L5_cnv_filt[1:4,i+4]))
}
t_result[,2] = p.adjust(t_result[,1], method = "fdr")
colnames(t_result) = c("p_value","FDR","log2FC")
rownames(t_result) = colnames(rate_L5_cnv_filt[,-(1:4)])
t_sig = t_result[t_result[,2] < 0.05,]

write.csv(t_sig, "Output_data/Canavanine_Taxon_L5_t_test_significant.csv")
write.csv(t_result, "Output_data/Canavanine_Taxon_L5_t_test_result.csv")

# Calculate Mean and SD of relative abundance of significantly different bacteria
rate_L5_cnv_sig = cbind(rate_L5_cnv_filt[,1:4],
                        rate_L5_cnv_filt[,c(rownames(t_sig[t_sig[,3] > 0,]), rownames(t_sig[t_sig[,3] < 0,]))])
rate_L5_cnv_sig_ave = rate_L5_cnv_sig[seq(1,16,4),1:4]
ave = rbind(colMeans(rate_L5_cnv_sig[1:4,-(1:4)]), colMeans(rate_L5_cnv_sig[5:8,-(1:4)]),
            colMeans(rate_L5_cnv_sig[9:12,-(1:4)]), colMeans(rate_L5_cnv_sig[13:16,-(1:4)]))
rate_L5_cnv_sig_ave = cbind(rate_L5_cnv_sig_ave, ave)

sd = rbind(apply(rate_L5_cnv_sig[1:4,-(1:4)],2,sd), apply(rate_L5_cnv_sig[5:8,-(1:4)],2,sd),
           apply(rate_L5_cnv_sig[9:12,-(1:4)],2,sd), apply(rate_L5_cnv_sig[13:16,-(1:4)],2,sd))
rate_L5_cnv_sig_sd = cbind(rate_L5_cnv_sig_ave[,1:4],sd)

rate_L5_cnv_sig_ave_2 = rate_L5_cnv_sig_ave[order(rate_L5_cnv_sig_ave$Conc),]
rate_L5_cnv_sig_sd_2 = rate_L5_cnv_sig_sd[order(rate_L5_cnv_sig_sd$Conc),]

# Prepare taxonomy name table
taxon_name = str_split(colnames(rate_L5_cnv_sig_ave_2)[-(1:4)], pattern = ";", simplify = T)
taxon_name = cbind(taxon_name,NA)
taxon_name[,6] = paste(str_sub(taxon_name[,3],start = 6),"(c); ",
                       str_sub(taxon_name[,4],start = 6),"(o); ",
                       str_sub(taxon_name[,5],start = 6),"(f)")
taxon_name[c(13:15,18),6] = paste(str_sub(taxon_name[c(13:15,18),2],start = 6),"(p); ",
                                  str_sub(taxon_name[c(13:15,18),3],start = 6),"(c); ",
                                  str_sub(taxon_name[c(13:15,18),4],start = 6),"(o)")

# Show relative abundance as barplot
barplot(rate_L5_cnv_sig_ave_2[,5])

file_name = "Output_data/Canavanine_Taxon_L5_significant_barplot.pdf"
pdf(file_name, width=12, height=10.0)
par(mfrow = c(4,6))
for (i in 1:(ncol(rate_L5_cnv_sig_ave_2)-4)) {
  xm = rate_L5_cnv_sig_ave_2[,i+4]
  xs = rate_L5_cnv_sig_sd_2[,i+4]
  bp = barplot(xm, names.arg = c("C","L","M","H"), col = adjustcolor(cols13[c(rep(5,12),rep(2,10))][i],alpha.f = 0.6), las = 1,
               ylim = c(0,max(xm + xs)*1.05), ylab = "Relative frequency (%)")
  arrows(bp,xm+xs,bp,xm-xs, angle = 90, code = 3, length = 0.01)
  box(lty = 1)
  mtext(side = 3, text = paste(str_split(taxon_name[i,6],";",simplify = T)[1],"\n",
                               str_split(taxon_name[i,6],";",simplify = T)[2],"\n",
                               str_split(taxon_name[i,6],";",simplify = T)[3]),
        line = 0.5, cex = 0.7)
}
plot(x=1, y=1, ann = F, xaxt="n", yaxt="n", bty = "n", ty = "n")
legend("left", legend = c("C: Control", "L: Low", "M: Middle", "H: High"), col = NA, box.col = NA)
dev.off()

# Show fold-change as heatmap
rate_L5_cnv_sig_ave_2 = cbind(rate_L5_cnv_sig_ave_2[,1:4], log2(rate_L5_cnv_sig_ave_2[,-(1:4)]/rate_L5_cnv_sig_ave_2[c(1,1,1,1),-(1:4)]))
heatmap(t(rate_L5_cnv_sig_ave_2[,-(1:4)]))

file_name = "Output_data/Canavanine_Taxon_L5_significant_Heatmap_FC.pdf"
pdf(file_name, width=8.0, height=6)
heatmap.2(t(rate_L5_cnv_sig_ave_2[,-(1:4)]),breaks = seq(-2,2,length.out = 101),
          key.title = "", keysize = 1.4, na.color="gray",dendrogram = "row",Colv = NA,Rowv = T,scale = "n",key.xlab = "log2 Fold change",
          density.info = "n", trace="n", col = colorRampPalette(c("blue", "white","red"))(100),srtCol = 0, cexRow = 0.8,
          labRow = taxon_name[,6], labCol = c("C","L","M","H"),
          margins = c(3,27), lwid = c(1.2,4), lhei = c(1,4),cexCol = 1,
          colsep = c(),rowsep = c(12),sepcolor = "black",sepwidth = c(0.01,0.01),main = "")
mtext(side = 3, text = "Differential Abundance Bacteria Family\n(Control vs High, FDR < 0.05)", cex = 1.2, font = 2, line = 0)
dev.off()


# Taxonomy, Species level
taxon_L7 = read.csv("Input_data/level-7_s132.csv", header = T, as.is = T)
taxon_L7_head = read.csv("Input_data/level-7_s132.csv", header = F, as.is = T)[1,]
colnames(taxon_L7) = taxon_L7_head
numcol = ncol(taxon_L7)
taxon_L7_cnv = cbind(taxon_L7[,c(1, (numcol-2):numcol)], 
                     taxon_L7[,2:(numcol-3)][,apply(taxon_L7[,2:(numcol-3)],2,max) > 0])

rate_L7_cnv = cbind(taxon_L7_cnv[,1:4], taxon_L7_cnv[,-(1:4)]/rowSums(taxon_L7_cnv[,-(1:4)])*100)

# I focused on Bacillaceae family, because Bacillaceae was the most abundant among significantly different bacteria, and
# Bacillus are known to possess "L-arginine degradation pathway" (picrust result, later).
# Extract relative abundance of Bacillaceae family at Species level
rate_L7_baci = cbind(rate_L7_cnv[,1:4],rate_L7_cnv[,grep("Bacillaceae",colnames(rate_L7_cnv))])
rate_L7_baci = cbind(rate_L7_baci[,1:4],rate_L7_baci[,-(1:4)][,colSums(rate_L7_baci[,-(1:4)])>0])

# Check taxonomy name and Sum relative abundance of unidentified species
# (e.g. "uncultured bacterium" + "uncultured compost bacterium" + "unidentified" --> "unidentified")
colnames(rate_L7_baci)
rate_L7_baci_2 = cbind(rate_L7_baci[,1:12], rowSums(rate_L7_baci[,13:16]), rowSums(rate_L7_baci[,17:18]),
                       rate_L7_baci[,19:25], rowSums(rate_L7_baci[,26:27]), rowSums(rate_L7_baci[,28:29]),
                       rate_L7_baci[,30:31], rowSums(rate_L7_baci[,32:35]), rowSums(rate_L7_baci[,36:40]))
colnames(rate_L7_baci_2) = colnames(rate_L7_baci)[c(1:12,16:17,19:25,27,29:31,35,40)]
colnames(rate_L7_baci_2) = gsub("D_6__uncultured bacterium","__",colnames(rate_L7_baci_2))
colnames(rate_L7_baci_2)

# Prepare taxonomy name table
baci_name = str_split(colnames(rate_L7_baci_2)[-(1:4)],";",simplify = T)
baci_name = cbind(baci_name,NA)
baci_name[c(2:8,11,13,15:17),8] = str_sub(baci_name[c(2:8,11,13,15:17),7], start = 6)
baci_name[-c(2:8,11,13,15:17),8] = paste(str_sub(baci_name[-c(2:8,11,13,15:17),6], start = 6),"(g)")
baci_name[23,8] = paste(baci_name[23,7],"(g)")

# Calculate Mean and SD of relative abundance of significantly different bacteria
rate_baci_ave = cbind(rate_L7_baci_2[c(1,9,13,5),1:4],
                      rbind(apply(rate_L7_baci_2[1:4,-(1:4)],2,mean), apply(rate_L7_baci_2[9:12,-(1:4)],2,mean),
                            apply(rate_L7_baci_2[13:16,-(1:4)],2,mean), apply(rate_L7_baci_2[5:8,-(1:4)],2,mean)))
rate_baci_sd = cbind(rate_L7_baci_2[c(1,9,13,5),1:4],
                     rbind(apply(rate_L7_baci_2[1:4,-(1:4)],2,sd), apply(rate_L7_baci_2[9:12,-(1:4)],2,sd),
                           apply(rate_L7_baci_2[13:16,-(1:4)],2,sd), apply(rate_L7_baci_2[5:8,-(1:4)],2,sd)))

# Show composition of Bacillaceae family at Species level as stacked barplot
barplot(t(rate_baci_ave[,rev(5:27)]))

col_22 = c(colorRampPalette(c(cols13[1],"white"))(4)[c(2)], colorRampPalette(c(cols13[2],"white"))(8)[c(7,5,3,1)],
           colorRampPalette(c(cols13[3],"white"))(8)[c(7,5,3,1)], colorRampPalette(c(cols13[1],"white"))(10)[c(9,7,5,3,1)],
           colorRampPalette(c(cols13[7],"white"))(8)[c(7,5,3,1)], colorRampPalette(c(cols13[10],"white"))(6)[c(5,3,1)],
           cols13[4])

file_name = "Output_data/Canavanine_Taxon_L7_Bacillaceae_barplot.pdf"
pdf(file_name, width=6.0, height= 4)
layout(matrix(c(1,1,1,2,2),nrow = 1))
bp = barplot(t(rate_baci_ave[,rev(5:27)]),ylab = "Relative Frequency (%)", las = 1, names = c("C","L","M","H"), main = "", 
             col = adjustcolor(rev(c(col_22,gcol[10])),alpha.f = 0.8), border = T, 
             ylim = c(0,max(rowSums(rate_baci_ave[1:4,rev(5:27)]))*1.05), xlab = "Canavanine-treatment")
box(bty = "l")
plot(x=1, y=1, ann = F, xaxt="n", yaxt="n", bty = "n", ty = "n")
legend("left", legend = baci_name[,8], fill = adjustcolor(c(col_22,gcol[10]), alpha.f = 0.8),
       border = T,cex = 0.9, xpd = T, title = "Bacillaceae (Family)")
dev.off()


# Abundance at ASV level
asv_table = read.table("Input_data/exported-table/asv_table.txt", sep = "\t", header = F, as.is = T)
taxonomy = read.table("Input_data/exported-taxonomy/taxonomy.tsv", sep = "\t", header = T, as.is = T)
colnames(asv_table) = c("ASV_ID",taxon_L7$index)
asv_cnv = asv_table

# Extract relative abundance of Bacillus genus at ASV level
asv_bacillus = subset(asv_cnv, ASV_ID %in% taxonomy[grep("D_5__Bacillus", taxonomy$Taxon),1])
asv_bacillus = asv_bacillus[rowSums(asv_bacillus[,-1]) > 0,]
rate_bacillus = t(asv_bacillus[,-1])/colSums(asv_cnv[,-1])*100
colnames(rate_bacillus) = asv_bacillus$ASV_ID

rate_bacillus = rate_bacillus[c(1:4,9:16,5:8),]
rate_bacillus_ave = cbind(rate_L7_baci_2[c(1,9,13,5),1:4],
                          rbind(apply(rate_bacillus[1:4,],2,mean), apply(rate_bacillus[5:8,],2,mean),
                                apply(rate_bacillus[9:12,],2,mean), apply(rate_bacillus[13:16,],2,mean)))

rate_bacillus_ave_2 = rate_bacillus_ave[,c(1:4,c(5:76)[order(colMeans(rate_bacillus_ave[,-(1:4)]),decreasing = T)])]

# Show composition of Bacillus genus at ASV level as stacked barplot
barplot(t(rate_bacillus_ave_2[,rev(5:76)]))

file_name = "Output_data/Canavanine_ASV_Bacillus_barplot.pdf"
pdf(file_name, width=6.0, height= 4)
layout(matrix(c(1,1,1,2,2),nrow = 1))
bp = barplot(t(rate_bacillus_ave_2[,rev(5:76)]),ylab = "Relative Frequency (%)", las = 1, names = c("C","L","M","H"), main = "", 
             col = adjustcolor(rev(c(cols13[c(1:7,10:11,13)],gcol[rep(10,62)])),alpha.f = 0.8), border = T, 
             ylim = c(0,max(rowSums(rate_bacillus_ave_2[1:4,rev(5:76)]))*1.05), xlab = "Canavanine-treatment")
box(bty = "l")
plot(x=1, y=1, ann = F, xaxt="n", yaxt="n", bty = "n", ty = "n")
legend("left", legend = paste("ASV",c(1:10,"11-72")), fill = adjustcolor(c(cols13[c(1:7,10:11,13)],gcol[10]), alpha.f = 0.8),
       border = T,cex = 0.9, xpd = T, title = "Bacillus (Genus)")
dev.off()


# Taxonomy, Phylum level
taxon_L2 = read.csv("Input_data/level-2_s132.csv", header = T, as.is = T)
taxon_L2_head = read.csv("Input_data/level-2_s132.csv", header = F, as.is = T)[1,]
colnames(taxon_L2) = taxon_L2_head
numcol = ncol(taxon_L2)
taxon_L2_cnv = cbind(taxon_L2[,c(1, (numcol-2):numcol)], 
                     taxon_L2[,2:(numcol-3)][,apply(taxon_L2[,2:(numcol-3)],2,max) > 0])

rate_L2_cnv = cbind(taxon_L2_cnv[,1:4], taxon_L2_cnv[,-(1:4)]/rowSums(taxon_L2_cnv[,-(1:4)])*100)

# Extract relative abundance of dominant phylum
rate_L2_cnv_rm1 = cbind(rate_L2_cnv[,1:4],rate_L2_cnv[,-(1:4)][,colMeans(rate_L2_cnv[,-(1:4)]) > 1], 
                        others = rowSums(rate_L2_cnv[,-(1:4)][,colMeans(rate_L2_cnv[,-(1:4)]) < 1]))

ave = rbind(colMeans(rate_L2_cnv_rm1[1:4,-(1:4)]), colMeans(rate_L2_cnv_rm1[5:8,-(1:4)]),
            colMeans(rate_L2_cnv_rm1[9:12,-(1:4)]), colMeans(rate_L2_cnv_rm1[13:16,-(1:4)]))
rate_L2_cnv_ave = cbind(rate_L2_cnv_rm1[seq(1,16,4),1:4],ave)
rate_L2_cnv_ave = rate_L2_cnv_ave[order(rate_L2_cnv_ave$Conc),]

# Show composition of dominant phyla as barplot
barplot(t(rate_L2_cnv_ave[,rev(5:15)]))

file_name = "Output_data/Canavanine_Taxon_L2_barplot.pdf"
pdf(file_name, width=8.0, height=5)
par(mfrow = c(1,2))
barplot(t(rate_L2_cnv_ave[,rev(5:15)]),names = c("C","L","M","H"),col = adjustcolor(rev(c(cols13[1:10],gcol[10])),alpha.f = 0.8),
        ylim = c(0,100), ylab = "Relative Frequency (%)", las = 1)
box()
plot(x=1, y=1, ann = F, xaxt="n", yaxt="n", bty = "n", ty = "n")
legend("right",legend = colnames(rate_L2_cnv_ave)[-(1:4)], fill = adjustcolor(c(cols13[1:10],gcol[10]),alpha.f = 0.8),
       title = "Dominant Phylum (Top10)",xpd = T, cex = 0.9)
dev.off()



## PICRUSt2
# Pathway level
pathway_desc = read.table("Input_data/description_mapfiles/metacyc_pathways_info.txt", sep = "\t",as.is = T,quote = "")
colnames(pathway_desc) = c("PATH_ID","Description")

path_abund_cnv = read.table("Input_data/pathabun_exported/feature-table.biom.tsv", sep = "\t",as.is = T)
colnames(path_abund_cnv) = c("PATH_ID",metadata_cnv$Sample_ID)
path_abund_cnv = merge(pathway_desc,path_abund_cnv,by = "PATH_ID")

# Differential abundant test, ALDEx2
count_cnv = cbind(path_abund_cnv[,1:2], round(path_abund_cnv[,-(1:2)]))
rownames(count_cnv) = count_cnv$PATH_ID
conds = c(rep("C",4),rep("H",4))
x_cnv = aldex(count_cnv[,c(3:6,7:10)], conds, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=FALSE)
cnv_sig = x_cnv[x_cnv$we.eBH < 0.05,]

# Transform count --> Relative abundance (%)
# Calculate Mean of relative abundance of pathways
path_rate_cnv = cbind(path_abund_cnv[,1:2], t(t(path_abund_cnv[,-(1:2)])/colSums(path_abund_cnv[,-(1:2)])*100))

path_rate_cnv_ave = cbind(path_rate_cnv[,1:2],matrix(NA,nrow = nrow(path_rate_cnv), ncol = 4))
for (i in 1:4) {
  path_rate_cnv_ave[,i+2] = rowMeans(path_rate_cnv[,(4*i-1):(4*i+2)])
}
colnames(path_rate_cnv_ave)[3:6] = colnames(path_rate_cnv)[seq(3,18,4)]
path_rate_cnv_ave = path_rate_cnv_ave[,c(1:3,5:6,4)]

# Extract relative abundance of significantly different pathways
path_rate_cnv_sig = subset(path_rate_cnv_ave,PATH_ID %in% rownames(cnv_sig))
path_rate_cnv_sig_fc = cbind(path_rate_cnv_sig[,1:2],log2(path_rate_cnv_sig[,3:6]/path_rate_cnv_sig[,3]))
path_rate_cnv_sig_fc$Description = gsub("&alpha;","alpha",path_rate_cnv_sig_fc$Description)
path_rate_cnv_sig_fc$Description = gsub("&beta;","beta",path_rate_cnv_sig_fc$Description)


# Show fold-change of significant pathways as heatmap
file_name = "Output_data/Picrust_Enriched_pathway_Heatmap.pdf"
pdf(file_name, width=9.0, height=9)
heatmap.2(as.matrix(path_rate_cnv_sig_fc[path_rate_cnv_sig_fc$`canavanine-high1`>0, 3:6]),
          key.xlab = "log2 Fold change", key.title = "", keysize = 1.4,na.color="black",dendrogram = "row",Colv = NA,
          density.info = "n", trace="n", col = colorRampPalette(c("blue", "white","red"))(100),srtCol = 0, 
          labRow = path_rate_cnv_sig_fc[path_rate_cnv_sig_fc$`canavanine-high1`>0,]$Description,
          labCol = c("C","L","M","H"),cexCol = 1,breaks = seq(-1.5,1.5,length.out = 101),
          margins = c(3,28), lwid = c(1.2,4),lhei = c(1,5.7),
          colsep = c(),rowsep = c(),sepcolor = "black",sepwidth = c(0.01,0.01),main = "")
mtext(side = 3, text = "Enriched pathway\n(Canavanine, Control vs High, FDR < 0.05)", cex = 1.5, font = 2, line = 0)
dev.off()

file_name = "Output_data/Picrust_Depleted_pathway_Heatmap.pdf"
pdf(file_name, width=9.0, height=9)
heatmap.2(as.matrix(path_rate_cnv_sig_fc[path_rate_cnv_sig_fc$`canavanine-high1`<0, 3:6]),
          key.xlab = "log2 Fold change", key.title = "", keysize = 1.4,na.color="black",dendrogram = "row",Colv = NA,
          density.info = "n", trace="n", col = colorRampPalette(c("blue", "white","red"))(100),srtCol = 0, 
          labRow = path_rate_cnv_sig_fc[path_rate_cnv_sig_fc$`canavanine-high1`<0,]$Description,
          labCol = c("C","L","M","H"),cexCol = 1,breaks = seq(-1.5,1.5,length.out = 101),
          margins = c(3,28), lwid = c(1.2,4),lhei = c(1,5.7),
          colsep = c(),rowsep = c(),sepcolor = "black",sepwidth = c(0.01,0.01),main = "")
mtext(side = 3, text = "Depleted pathway\n(Canavanine, Control vs High, FDR < 0.05)", cex = 1.5, font = 2, line = 0)
dev.off()

# Output list of significantly different pathways
write.csv(subset(path_rate_cnv_sig_fc,PATH_ID %in% rownames(cnv_sig[cnv_sig$diff.btw>0,]))[,1:2],"Output_data/Picrust_Enriched_Pathway_list.csv",row.names = F)
write.csv(subset(path_rate_cnv_sig_fc,PATH_ID %in% rownames(cnv_sig[cnv_sig$diff.btw<0,]))[,1:2],"Output_data/Picrust_Depleted_Pathway_list.csv",row.names = F)


## Picrust2
# Enzyme level
ec_info = read.table("Input_data/description_mapfiles/ec_level4_info.tsv", sep = "\t",as.is = T, quote = "")
ec_abund_cnv = read.table("Input_data/ec_exported/feature-table.biom.tsv", sep = "\t",as.is = T)
colnames(ec_abund_cnv) = c("EC",metadata_cnv$Sample_ID)
ec_abund_cnv = merge(ec_info,ec_abund_cnv, by.x = "V1", by.y = "EC")
colnames(ec_abund_cnv)[1:2] = c("EC","EC_desc")

# Differential abundant test, ALDEx2
ec_count_cnv = cbind(ec_abund_cnv[,1:2], round(ec_abund_cnv[,-(1:2)]))
rownames(ec_count_cnv) = ec_count_cnv$EC
conds = c(rep("C",4),rep("H",4))
x_ec_cnv = aldex(ec_count_cnv[,c(3:6,7:10)], conds, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=FALSE)
ec_cnv_sig = x_ec_cnv[x_ec_cnv$we.eBH < 0.05,]

# Output list of significantly different pathways
write.csv(subset(ec_abund_cnv,EC %in% rownames(ec_cnv_sig[ec_cnv_sig$diff.btw>0,]))[,1:2],"Output_data/Picrust_Enriched_EC_list.csv",row.names = F)
write.csv(subset(ec_abund_cnv,EC %in% rownames(ec_cnv_sig[ec_cnv_sig$diff.btw<0,]))[,1:2],"Output_data/Picrust_Depleted_EC_list.csv",row.names = F)


## Picrust2
# KEGG Ortholog level
ko_info = read.table("Input_data/description_mapfiles/ko_info.tsv", sep = "\t",as.is = T, quote = "")
ko_abund_cnv = read.table("Input_data/ko_exported/feature-table.biom.tsv", sep = "\t",as.is = T)
colnames(ko_abund_cnv) = c("KO",metadata_cnv$Sample_ID)
ko_abund_cnv = merge(ko_info,ko_abund_cnv, by.x = "V1", by.y = "KO")
colnames(ko_abund_cnv)[1:2] = c("KO","KO_desc")

# Differential abundant test, ALDEx2
ko_count_cnv = cbind(ko_abund_cnv[,1:2], round(ko_abund_cnv[,-(1:2)]))
rownames(ko_count_cnv) = ko_count_cnv$KO
conds = c(rep("C",4),rep("H",4))
x_ko_cnv = aldex(ko_count_cnv[,c(3:6,7:10)], conds, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=FALSE)
ko_cnv_sig = x_ko_cnv[x_ko_cnv$we.eBH < 0.05,]

# Output list of significantly different pathways
write.csv(subset(ko_abund_cnv,KO %in% rownames(ko_cnv_sig[ko_cnv_sig$diff.btw>0,]))[,1:2],"Output_data/Picrust_Enriched_KEGG_Orthology_list.csv",row.names = F)
write.csv(subset(ko_abund_cnv,KO %in% rownames(ko_cnv_sig[ko_cnv_sig$diff.btw<0,]))[,1:2],"Output_data/Picrust_Depleted_KEGG_Orthology_list.csv",row.names = F)


