#install.packages("BiocManager")
#BiocManager::install("coRdon")
#install.packages("BiocManager")
#install.packages("rphast")
library(rphast)
library(coRdon)

# formatting ###################################################################################################################################################

aln<-read.msa('OG0002887.ref.mafft.n.aln')
tre<-read.newick.tree("sp.nwk")
#feats<- read.feat("aln/OG1.cds.fas.transdecoder_dir/longest_orfs.cds.best_candidates.gff3")
#cds<-extract.feature.msa(aln, feats[feats$feature=="CDS",], do4d = FALSE, pointer.only = FALSE)
#fds<-extract.feature.msa(aln, do4d = TRUE, pointer.only = FALSE)
#utr<-extract.feature.msa(aln, feats[feats$feature=="UTR",], do4d = FALSE, pointer.only = FALSE)

# substitution mapping in cds ##################################################################################################################################

model<-phyloFit(aln, tree=tre, subst.mod="REV")
BAT.result<-classify.muts.bgc (aln,model,"BAT")
BRO.result<-classify.muts.bgc (aln,model,"BRO")
BGM.result<-classify.muts.bgc (aln,model,"BGM")
#write.table(dmag.result,file="out.tsv",quote=FALSE, sep ="\t")

# bgGC model fitting ###########################################################################################################################################

bgc.BAT.init.model<-add.ls.mod(model,"BAT",separate.params="bgc[0,2000]")
bgc.BAT.model<-phyloFit(aln,init.mod=bgc.BAT.init.model)
LRT <- -2*(bgc.BAT.init.model$likelihood-bgc.BAT.model$likelihood)
p.value <- 1-pchisq(LRT,df=1)
if (p.value < 0.001) {
  bgc.BAT <- bgc.BAT.model$ls.model$bgc
  } else {
  bgc.BAT <- NA
  }
bgc.BAT

bgc.BRO.init.model<-add.ls.mod(model,"BRO",separate.params="bgc[0,2000]")
bgc.BRO.model<-phyloFit(aln,init.mod=bgc.BRO.init.model)
LRT <- -2*(bgc.BRO.init.model$likelihood-bgc.BRO.model$likelihood)
p.value <- 1-pchisq(LRT,df=1)
if (p.value < 0.001) {
  bgc.BRO <- bgc.BRO.model$ls.model$bgc
} else {
  bgc.BRO <- NA
}
bgc.BRO

bgc.BGM.init.model<-add.ls.mod(model,"BGM",separate.params="bgc[0,2000]")
bgc.BGM.model<-phyloFit(aln,init.mod=bgc.BGM.init.model)
LRT <- -2*(bgc.BGM.init.model$likelihood-bgc.BGM.model$likelihood)
p.value <- 1-pchisq(LRT,df=1)
if (p.value < 0.001) {
  bgc.BGM <- bgc.BGM.model$ls.model$bgc
} else {
  bgc.BGM <- NA
}
bgc.BGM

# GC of cds and fds  ###########################################################################################################################################

gc.cds.BAT<-gc.content.msa(aln, seq = "BAT", ignore.missing = TRUE, ignore.gaps = TRUE)
gc.cds.BRO<-gc.content.msa(aln, seq = "BRO", ignore.missing = TRUE, ignore.gaps = TRUE)
gc.cds.BGM<-gc.content.msa(aln, seq = "BGM", ignore.missing = TRUE, ignore.gaps = TRUE)

#gc.fds.BAT<-gc.content.msa(fds, seq = "BAT", ignore.missing = TRUE, ignore.gaps = TRUE)
#gc.fds.BRO<-gc.content.msa(fds, seq = "BRO", ignore.missing = TRUE, ignore.gaps = TRUE)
#gc.fds.BGM<-gc.content.msa(fds, seq = "BGM", ignore.missing = TRUE, ignore.gaps = TRUE)

#gc.utr.dmag<-gc.content.msa(utr, seq = "dmag", ignore.missing = TRUE, ignore.gaps = TRUE)
#gc.utr.dpul<-gc.content.msa(utr, seq = "dpul", ignore.missing = TRUE, ignore.gaps = TRUE)
#gc.utr.tcan<-gc.content.msa(utr, seq = "tcan", ignore.missing = TRUE, ignore.gaps = TRUE)


# pairwise distances of spp  ###################################################################################################################################

#tcan.tusa.dist<-pairwise.diff.msa(alignment, seq1 = "tcan", seq2 = "tusa", ignore.missing = TRUE,ignore.gaps = TRUE)
#dmag.dpul.dist<-pairwise.diff.msa(alignment, seq1 = "dmag", seq2 = "dpul", ignore.missing = TRUE,ignore.gaps = TRUE)
#lart.lubb.dist<-pairwise.diff.msa(alignment, seq1 = "lart", seq2 = "lubb", ignore.missing = TRUE,ignore.gaps = TRUE)

# codon usage bias  ############################################################################################################################################

aln <- readSet(file="OG0002887.ref.mafft.n.aln")
alnct <- codonTable(aln)

# enc ##########################################################################################################################################################

enc <- ENC(alnct)
enc.BAT<-enc[1]
enc.BRO<-enc[2]
enc.BGM<-enc[3]

# scuo #########################################################################################################################################################

scuo <- SCUO(alnct)
scuo.BAT<-enc[1]
scuo.BRO<-enc[2]
scuo.BGM<-enc[3]

