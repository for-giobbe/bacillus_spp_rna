echo "OG BAT_W-S BAT_S-W BAT_W-W BAT_S-S BRO_W-S BRO_S-W  BRO_W-W BRO_S-S BGM_W-S BGM_S-W BGM_W-W BGM_S-S BGC_BAT BGC_BAT_pval BGC_BRO BGC_BRO_pval BGC_BGM BGC_BGM_pval ENC_BAT ENC_BRO ENC_BGM SCUO_BAT SCUO_BRO SCUO_BGM MILC_BAT MILC_BRO MILC_BGM BAT_4FDS BRO_4FDS BGM_4FDS" > substiution_mapping.tab

printf '#install.packages("BiocManager") \n' > r_blueprint.tmp
printf '#BiocManager::install("coRdon") \n' >> r_blueprint.tmp
printf '#install.packages("BiocManager") \n' >> r_blueprint.tmp
printf '#install.packages("rphast") \n' >> r_blueprint.tmp
printf 'library(rphast) \n' >> r_blueprint.tmp
printf 'library(coRdon) \n' >> r_blueprint.tmp

printf '# formatting ################################################################################################################################################### \n' >> r_blueprint.tmp

printf 'aln<-read.msa("substitute_this") \n' >> r_blueprint.tmp
printf 'tre<-read.newick.tree("sp.nwk") \n' >> r_blueprint.tmp
printf 'feats<-read.feat("gff3.tmp") \n' >> r_blueprint.tmp
printf 'feats \n' >> r_blueprint.tmp
printf 'fds<-extract.feature.msa(aln,features=feats,do4d=TRUE,pointer.only=FALSE) \n' >> r_blueprint.tmp
printf 'fds \n' >> r_blueprint.tmp
printf '# substitution mapping in cds ################################################################################################################################## \n' >> r_blueprint.tmp

printf 'model<-phyloFit(aln, tree=tre, subst.mod="REV") \n' >> r_blueprint.tmp
printf 'BAT.result<-classify.muts.bgc (aln,model,"BAT") \n' >> r_blueprint.tmp
printf 'BRO.result<-classify.muts.bgc (aln,model,"BRO") \n' >> r_blueprint.tmp
printf 'BGM.result<-classify.muts.bgc (aln,model,"BGM") \n' >> r_blueprint.tmp

printf '# bgGC model fitting ########################################################################################################################################### \n' >> r_blueprint.tmp

printf 'bgc.BAT.init.model<-add.ls.mod(model,"BAT",separate.params="bgc[0,2000]") \n' >> r_blueprint.tmp
printf 'bgc.BAT.model<-phyloFit(aln,init.mod=bgc.BAT.init.model) \n' >> r_blueprint.tmp
printf 'LRT <- -2*(bgc.BAT.init.model$likelihood-bgc.BAT.model$likelihood) \n' >> r_blueprint.tmp
printf 'BAT.p.value <- 1-pchisq(LRT,df=1) \n' >> r_blueprint.tmp
printf 'if (BAT.p.value < 0.001) { \n' >> r_blueprint.tmp
printf '  bgc.BAT <- bgc.BAT.model$ls.model$bgc \n' >> r_blueprint.tmp
printf '  } else { \n' >> r_blueprint.tmp
printf '  bgc.BAT <- NA \n' >> r_blueprint.tmp
printf '  } \n' >> r_blueprint.tmp

printf 'bgc.BRO.init.model<-add.ls.mod(model,"BRO",separate.params="bgc[0,2000]") \n' >> r_blueprint.tmp
printf 'bgc.BRO.model<-phyloFit(aln,init.mod=bgc.BRO.init.model) \n' >> r_blueprint.tmp
printf 'LRT <- -2*(bgc.BRO.init.model$likelihood-bgc.BRO.model$likelihood) \n' >> r_blueprint.tmp
printf 'BRO.p.value <- 1-pchisq(LRT,df=1) \n' >> r_blueprint.tmp
printf 'if (BRO.p.value < 0.001) { \n' >> r_blueprint.tmp
printf '  bgc.BRO <- bgc.BRO.model$ls.model$bgc \n' >> r_blueprint.tmp
printf '} else { \n' >> r_blueprint.tmp
printf '  bgc.BRO <- NA \n' >> r_blueprint.tmp
printf '} \n' >> r_blueprint.tmp

printf 'bgc.BGM.init.model<-add.ls.mod(model,"BGM",separate.params="bgc[0,2000]") \n' >> r_blueprint.tmp
printf 'bgc.BGM.model<-phyloFit(aln,init.mod=bgc.BGM.init.model) \n' >> r_blueprint.tmp
printf 'LRT <- -2*(bgc.BGM.init.model$likelihood-bgc.BGM.model$likelihood) \n' >> r_blueprint.tmp
printf 'BGM.p.value <- 1-pchisq(LRT,df=1) \n' >> r_blueprint.tmp
printf 'if (BGM.p.value < 0.001) { \n' >> r_blueprint.tmp
printf '  bgc.BGM <- bgc.BGM.model$ls.model$bgc \n' >> r_blueprint.tmp
printf '} else { \n' >> r_blueprint.tmp
printf '  bgc.BGM <- NA \n' >> r_blueprint.tmp
printf '} \n' >> r_blueprint.tmp

printf '# GC of cds and fds  ########################################################################################################################################### \n' >> r_blueprint.tmp

printf 'gc.fds.BAT<-gc.content.msa(fds, seq = "BAT", ignore.missing = TRUE, ignore.gaps = TRUE) \n' >> r_blueprint.tmp
printf 'gc.fds.BRO<-gc.content.msa(fds, seq = "BRO", ignore.missing = TRUE, ignore.gaps = TRUE) \n' >> r_blueprint.tmp
printf 'gc.fds.BGM<-gc.content.msa(fds, seq = "BGM", ignore.missing = TRUE, ignore.gaps = TRUE) \n' >> r_blueprint.tmp

printf '# codon usage bias  ############################################################################################################################################ \n' >> r_blueprint.tmp

printf 'aln <- readSet(file="substitute_this") \n' >> r_blueprint.tmp
printf 'alnct <- codonTable(aln) \n' >> r_blueprint.tmp

printf '# enc ########################################################################################################################################################## \n' >> r_blueprint.tmp

printf 'enc <- ENC(alnct, filtering = "hard",len.threshold = 80) \n' >> r_blueprint.tmp
printf 'enc.BAT<-enc[1] \n' >> r_blueprint.tmp
printf 'enc.BRO<-enc[2] \n' >> r_blueprint.tmp
printf 'enc.BGM<-enc[3] \n' >> r_blueprint.tmp

printf '# scuo ######################################################################################################################################################### \n' >> r_blueprint.tmp

printf 'scuo <- SCUO(alnct, filtering = "hard",len.threshold = 80) \n' >> r_blueprint.tmp
printf 'scuo.BAT<-scuo[1] \n' >> r_blueprint.tmp
printf 'scuo.BRO<-scuo[2] \n' >> r_blueprint.tmp
printf 'scuo.BGM<-scuo[3] \n' >> r_blueprint.tmp

printf '# milc ######################################################################################################################################################### \n' >> r_blueprint.tmp

printf 'milc <- MILC(alnct, ribosomal = FALSE, filtering = "hard", len.threshold = 80) \n' >> r_blueprint.tmp
printf 'milc.BAT<-milc[1] \n' >> r_blueprint.tmp
printf 'milc.BRO<-milc[2] \n' >> r_blueprint.tmp
printf 'milc.BGM<-milc[3] \n' >> r_blueprint.tmp

printf '# write results ################################################################################################################################################ \n' >> r_blueprint.tmp

printf 'sink("result.tmp") \n' >> r_blueprint.tmp
printf 'cat(c(BAT.result$W.to.S, BAT.result$S.to.W, BAT.result$W.to.W, BAT.result$S.to.S, BRO.result$W.to.S, BRO.result$S.to.W, BRO.result$W.to.W, BRO.result$S.to.S, BGM.result$W.to.S, BGM.result$S.to.W, BGM.result$W.to.W, BGM.result$S.to.S, bgc.BAT, BAT.p.value, bgc.BRO, BRO.p.value, bgc.BGM, BGM.p.value, enc.BAT, enc.BRO, enc.BGM, scuo.BAT, scuo.BRO, scuo.BGM, milc.BAT, milc.BRO, milc.BGM,gc.fds.BAT,gc.fds.BRO,gc.fds.BGM)) \n' >> r_blueprint.tmp
printf 'cat("\n") \n' >> r_blueprint.tmp
printf 'sink() \n' >> r_blueprint.tmp

echo -e "BAT\ttrinity\tgene\t1\tlength\nBAT\ttrinity\tCDS\t1\tlength" > gff3_blueprint.tmp

echo -e "(((BRO,BGM),BAT),PHY);" > sp.nwk

for i in *aln;

	do

	OG=$( echo "$i " | awk -F "." '{print $1}')
	echo -n "$OG " >> substiution_mapping.tab
	sed "s/substitute_this/$i/g" r_blueprint.tmp > r.tmp;

	length=$(expr -1 + $(grep -A 1 "BAT" $i | tail -1 | wc -c))
	sed "s/length/$length/g" gff3_blueprint.tmp > gff3.tmp
	Rscript r.tmp;
	cat result.tmp >> substiution_mapping.tab
	done

rm *tmp
