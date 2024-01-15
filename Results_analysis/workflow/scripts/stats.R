library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(UpSetR)
library(ggh4x)
library(stringr)

prefix <- snakemake@params[['prefix']]
outputpath <- snakemake@params[['outputpath']]
runname <- snakemake@params[['runname']]
threshold_avgmcc <- snakemake@params[['avgmcc']]
threshold_stdmcc <- snakemake@params[['stdmcc']]
classifnames <- snakemake@params[['classifnames']]

redundantfeat <- snakemake@input[['redundantfeat']]
infogainorder <- snakemake@input[['infogainorder']]
resultsstand <- snakemake@input[['resultsstand']]

setwd(outputpath)
dir.create(file.path(outputpath, 'STATS'), showWarnings = FALSE)

results_by_samp <- list()

for (run in prefix){
    runresult <- grep(paste(run,"_", sep=""), resultsstand, value = TRUE)
    runresult <- lapply(runresult, read.csv)
    runresult <- bind_rows(runresult)

    runresult$samp <- run

    results_by_samp[[run]] <- runresult
}

fullresults <- bind_rows(results_by_samp)
fullresults$typeofclassif <- sub("\\..*", "", fullresults$classifier)
fullresults <- fullresults %>% 
  rename(
    MCC = AVG_MCC
  )

for (classif in classifnames){
    fullresults$classifier <- gsub(classif$weka, classif$change, fullresults$classifier )
}
fullresults$classifier_strat <- paste(fullresults$classifier,fullresults$samp,sep="_")
fullresults$typeofclassif <- str_to_title(fullresults$typeofclassif)


fullresults.filtered <- subset(fullresults, MCC >= threshold_avgmcc & STD_MCC <= threshold_stdmcc)

write.csv(table(fullresults$classifier,fullresults$samp),'STATS/stat_classif_before_filters.csv')
write.csv(table(fullresults.filtered$classifier,fullresults.filtered$samp),paste('STATS/stat_classif_after_filters_MCC',gsub('\\.', '', as.character(threshold_avgmcc)),'_STD',gsub('\\.', '', as.character(threshold_stdmcc)),'.csv',sep=""))

## BOXPLOTS
order_classif <- fullresults[c("classifier","typeofclassif")]
order_classif <- order_classif[order(order_classif$typeofclassif),]
order_classif <- unique(order_classif$classif)


fullresults$classifier <- factor(x = fullresults$classifier, levels=order_classif )
supp.labs <- unique(fullresults$typeofclassif)

fullresults <- fullresults[!is.na(fullresults$MCC),]
fullresults <- fullresults[!is.na(fullresults$STD_MCC),]


fullresults$samp <- factor(fullresults$samp, 
                       levels = prefix )

boxplot.MCC<- ggplot(fullresults, aes(x=samp, y=MCC, color=classifier)) + 
  geom_boxplot()+
  facet_nested(.~typeofclassif+classifier,
  labeller = label_wrap_gen(5),
  strip = strip_nested(size = "variable"),
  scales = "free",
  space='free'
  )+
  guides(fill = guide_legend(override.aes = list(size = 0.1)))+
  labs(color="Type of classifier")+
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits = c(0, 1))+
  xlab("Models ([algorithm]_[Sampling identifier])")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,face="plain"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        axis.title=element_text(size=18,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key = element_rect(fill = NA),
        strip.text.x = element_text(
          size = 15, color = "black"),
        strip.background = element_rect(
          color="black", size=1.5, linetype="solid"))+
  scale_color_manual(values=c("#af0648", "#d20811","#fc3b5d", "#f1c232","#04a424","#666cf0","#789cfb", "#78bdfb"),drop=FALSE)+
  geom_hline(aes(yintercept = 0.7),color="red", linetype="dashed")

boxplot.STDMCC<- ggplot(fullresults, aes(x=samp, y=STD_MCC, color=classifier)) + 
  geom_boxplot() +
  facet_nested(.~typeofclassif+classifier,
  labeller = label_wrap_gen(5),
  strip = strip_nested(size = "variable"),
  scales = "free",
  space='free'
  )+
  guides(colour = guide_legend(nrow = 2,title="Classifier"))+
  scale_y_continuous(breaks=seq(0,0.7,by=0.1))+
  xlab("Run")+
  ylab("STD MCC")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_text(size=12,face="plain"),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=10,face="plain"),
        axis.title=element_text(size=18,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12, margin = margin(r = 1.5, unit = 'cm'), hjust = 0),
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,1,1), "pt"))+
  scale_color_manual(values=c("#af0648", "#d20811","#fc3b5d", "#f1c232","#04a424","#666cf0","#789cfb", "#78bdfb"),drop=FALSE)+
  geom_hline(aes(yintercept = 0.1),color="red", linetype="dashed")

## save boxplot plot MMC + STD MCC
boxplots.MCC_and_STDMCC <- cowplot::plot_grid(boxplot.MCC, boxplot.STDMCC, 
                                             ncol = 1, rel_heights = c(0.6, 0.7),
                                             align = 'v', axis = 'lr')
ggsave("STATS/boxplot_MCC_and_STDMCC_params_sensi_onlyMCCopti_FB_DASH_LINE.png",boxplots.MCC_and_STDMCC, width = 12,height = 8)

## COUNT NBFEAT

#### geom count NB FEAT, on FILTERED MODELS ####
order_classif <- fullresults.filtered[c("classifier","typeofclassif")]
order_classif <- order_classif[order(order_classif$typeofclassif),]
order_classif <- unique(order_classif$classif)

fullresults.filtered$classifier <- factor(x = fullresults.filtered$classifier, levels=order_classif )
supp.labs <- unique(fullresults.filtered$typeofclassif)

my.aspect.ratio = length(unique(fullresults.filtered$nbrOfFeatures))

fullresults.filtered$samp <- factor(fullresults.filtered$samp, 
                       levels = prefix )

geomcountplot.full.nbfeat <- ggplot(fullresults.filtered, aes(x=samp, y=nbrOfFeatures,fill=classifier)) +
  geom_count(aes(size = after_stat(prop), group = samp),stroke = 0.2,shape = 21)+
  facet_nested(.~typeofclassif+classifier,
  labeller = label_wrap_gen(5),
  strip = strip_nested(size = "variable"),
  scales = "free",
  space='free'
  )+
  scale_size_area(max_size = 4)+
  guides(fill = guide_legend(nrow = 2,override.aes = list(size = 5)))+
  labs(fill="Classifier", size="Proportion by classifier")+
  ylab("Signature length (# features)")+
  xlab("Run")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_text(size=10,face="plain"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 10, face = "plain"),
        axis.title=element_text(size=14,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(-0.1, "cm"),
        legend.box="vertical",
        legend.position = "bottom",
        strip.text.x = element_text(
          size = 12, color = "black"),
        strip.background = element_rect(
          color="black", size=1, linetype="solid"),
          legend.key = element_rect(fill = NA))+
  scale_y_continuous(breaks=seq(0,max(fullresults.filtered$nbrOfFeatures),1),limits=c(0,max(fullresults.filtered$nbrOfFeatures)),expand = c(0,0))+
  scale_fill_manual(values=c("#af0648", "#d20811","#fc3b5d", "#f1c232","#04a424","#666cf0","#789cfb", "#78bdfb"),drop=FALSE)
ggsave(paste('STATS/MCC',gsub('\\.', '', as.character(threshold_avgmcc)),'_STD',gsub('\\.', '', as.character(threshold_stdmcc)),'FILTERED_geomcountplot_proportionate_nbfeat.png',sep=""),geomcountplot.full.nbfeat,width = 12,height = 7)


## Stats mean/std/IQR

#### STATS perf (AVG_MCC, STD_MCC) AND signature length - before AND after filtering ####
stats.perfandsign <- setDT(fullresults)[,list(
                                  AVG_AVG_MCC=mean(MCC),
                                  STD_AVG_MCC=sd(MCC),
                                  IQR_AVG_MCC=IQR(MCC),
                                  AVG_STD_MCC=mean(STD_MCC),
                                  STD_STD_MCC=sd(STD_MCC),
                                  IQR_STD_MCC=IQR(STD_MCC),
                                  AVG_SIG_LENGTH=mean(nbrOfFeatures),
                                  STD_SIG_LENGTH=sd(nbrOfFeatures),
                                  IQR_SIG_LENGTH=IQR(nbrOfFeatures)),
                                 by=list(samp,classifier)]

stats.perfandsign.filtered <- setDT(fullresults.filtered)[,list(
  AVG_AVG_MCC=mean(MCC),
  STD_AVG_MCC=sd(MCC),
  IQR_AVG_MCC=IQR(MCC),
  AVG_STD_MCC=mean(STD_MCC),
  STD_STD_MCC=sd(STD_MCC),
  IQR_STD_MCC=IQR(STD_MCC),
  AVG_SIG_LENGTH=mean(nbrOfFeatures),
  STD_SIG_LENGTH=sd(nbrOfFeatures),
  IQR_SIG_LENGTH=IQR(nbrOfFeatures)),
  by=list(samp,classifier)]

write.csv(stats.perfandsign,'STATS/stats_perf_and_siglength_methods_nothreshold_ALLsamp.csv',row.names = TRUE)
write.csv(stats.perfandsign.filtered,paste('STATS/stats_perf_and_siglength_methods_ALLsamp_MCC',gsub('\\.', '', as.character(threshold_avgmcc)),'_STD',gsub('\\.', '', as.character(threshold_stdmcc)),'.csv',sep=""),row.names = TRUE)

## OVERLAPS FEATURES
##### open infogain ##### 
format_redundant_feats <- function(redundantfeat,infogainshort){
  redundantfeat$correlated_feature <- gsub("Transcripto_Gene__", "\\1", redundantfeat$correlated_feature)
  redundantfeat$target_feature <- gsub("Transcripto_Gene__", "\\1", redundantfeat$target_feature)
  temp <- colnames(infogainshort)
  infogainorder <- 1:length(temp)
  names(infogainorder) <- temp

  redundantfeat <- select(redundantfeat,c('target_feature','correlated_feature'))

  temp2 <- redundantfeat
  colnames(temp2) <- c('correlated_feature','target_feature')
  redundantfeat.full <- rbind(redundantfeat, temp2)

  redundantfeat.full <- redundantfeat.full %>% 
    group_by(target_feature) %>% 
    mutate(fullcorr = paste(correlated_feature, collapse= ','))

  redundantfeat.full$correlated_feature <- NULL
  redundantfeat.full <- redundantfeat.full[!duplicated(redundantfeat.full), ]
  redundantfeat.full$fullfeat <- apply( redundantfeat.full[ , c('target_feature','fullcorr') ] , 1 , paste , collapse = "," )
  redundantfeat.full$fullcorr <- NULL
  
  returnlist <- list('redundant'=redundantfeat.full,'infogainorder'=infogainorder)
}


#### check overlap in ONE sampling for models of ONE method ####
overlapmodels <- function(results,choicesamp,choicemethod,redundantfeat,infogainfeat){
  results$samp[]
  results <- subset(results, results$samp == choicesamp & results$classifier == choicemethod,select = c('ID','AttributeList'))
  nbmodels <- length(results$ID)
  results <- separate_rows(results, AttributeList,sep=',', convert = FALSE)
  results <- subset(results, ! results$AttributeList %in% c(1,length(infogainfeat))) 

  results$AttributeList <- names(infogainfeat)[match(results$AttributeList, infogainfeat)]
  results$AttributeList <- gsub("Transcripto_Gene__", "\\1",results$AttributeList)
  results <- merge(results,redundantfeat,by.x = 'AttributeList',by.y = 'target_feature',all.x=T)
  results <- results %>% 
    mutate(fullfeat = coalesce(fullfeat,AttributeList))
  results$AttributeList <- NULL
  results <- separate_rows(results, fullfeat,sep=',', convert = FALSE)
  results <- results[!duplicated(results), ]
  
  if(dim(results)[1] == 0) { return(results)}
  results <- setattr(results, "class", c("data.table"))

  results <- dcast(results, ID ~ fullfeat, fun.aggregate = function(x) 1L, fill = 0L)

  results.forupset <- results
  results <- t(results)
  colnames(results) <- results[1,]
  results <- results[-1,]
  results <- as.data.frame(results)

  results[] <- as.data.frame(sapply(results, as.numeric))
  results <- rowSums(results)
  results <- as.data.frame(results)
  results$feat <- row.names(results)
  results$percentmodels <- results$results*100/nbmodels
  results
  }

### get features + redundant features for filtered models by sampling
### INTRA-METHOD INTRA-SAMPLING
dir.create(file.path(outputpath, 'STATS/intramethod_intrasampling'), showWarnings = FALSE)

feat_by_samp_by_method <- list()
feat_by_method_by_samp <- list()
for (run in prefix){
    run_redundantfeat <- grep(paste("/",run,"/",sep=""), redundantfeat, value = TRUE)
    run_redundantfeat<- grep("FILTERED", run_redundantfeat, value = TRUE)
    run_redundantfeat <- read.csv(run_redundantfeat)

    run.igorder <- grep(paste("/",run,"/",sep=""), infogainorder, value = TRUE)
    run.igorder <- read.csv(run.igorder)

    runfeat <- format_redundant_feats(run_redundantfeat,run.igorder)

    for (classif in order_classif){
        for (searchedclassif in names(classifnames)){
            if (classifnames[[searchedclassif]]$change == classif){
                initialname <- searchedclassif

            }
        }
      run.modelsinmethod.overlaps <- overlapmodels(fullresults.filtered,run,classif,data.frame(runfeat[[1]]),runfeat[[2]])
      write_csv(run.modelsinmethod.overlaps,paste('STATS/intramethod_intrasampling/',classif,'_',run,'_FILTERED_MCC',threshold_avgmcc,'_STDMCC',threshold_stdmcc,'_involved_feats_overlaps.csv',sep=""))
      
      feat_by_samp_by_method[[run]][[classif]] <- run.modelsinmethod.overlaps$feat

      feat_by_method_by_samp[[classif]][[run]] <- run.modelsinmethod.overlaps$feat
    }
}

fromList <- function (input) {
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  row.names(data) <- elements
  return(data)
}

### UNI-METHOD INTER-SAMPLING
dir.create(file.path(outputpath, 'STATS/unimethod_intersampling'), showWarnings = FALSE)
for (method in names(feat_by_method_by_samp)){
  methodmatrix <- fromList(feat_by_method_by_samp[[method]])
  methodmatrix$intersect <- rowSums(methodmatrix) 
  write.csv(methodmatrix,paste('STATS/unimethod_intersampling/',method,'_MATRIX_unimethod_intersampling_feats_overlaps.csv',sep=""),row.names = TRUE)

  upsetPlot <- upset(fromList(feat_by_method_by_samp[[method]]), nintersects = NA)
  ggsave(paste('STATS/unimethod_intersampling/',method,'_UPSETPLOT_unimethod_intersampling_feats_overlaps.png',sep=""),ggplotify::as.ggplot(upsetPlot))
  
}

### INTER-METHOD INTRA-SAMPLING
dir.create(file.path(outputpath, 'STATS/intermethod_intrasampling'), showWarnings = FALSE)
maxoverlap_of_methods_by_sampling <- list()
for (samp in names(feat_by_samp_by_method)){
  sampmatrix <- fromList(feat_by_samp_by_method[[samp]])
  sampmatrix$intersect <- rowSums(sampmatrix) 
  write.csv(sampmatrix,paste('STATS/intermethod_intrasampling/',samp,'_MATRIX_intermethod_intrasampling_feats_overlaps.csv',sep=""),row.names = TRUE)
  upsetPlot <- upset(fromList(feat_by_samp_by_method[[samp]]),
                      nsets = length(feat_by_samp_by_method[[samp]]),
                      sets = rev(names(feat_by_samp_by_method[[samp]])),
                      keep.order = T,
                      nintersects = NA,
                      point.size = 8,
                      line.size = 2,
                      mb.ratio = c(0.48, 0.52),
                      matrix.dot.alpha = 0.6,
                      text.scale = 3.0)
  ggsave(paste('STATS/intermethod_intrasampling/',samp,'_UPSETPLOT_intermethod_intrasampling_feats_overlaps.png',sep=""),ggplotify::as.ggplot(upsetPlot),width = 20,height = 7)
  sampmatrix_for_interM_interS <- subset(sampmatrix['intersect'], sampmatrix$intersect == max(sampmatrix$intersect)) 

  sampmatrix_for_interM_interS <- row.names(sampmatrix_for_interM_interS)
  sampname <- paste('intersect_',samp,'_maxis',as.character(max(sampmatrix$intersect)),sep="")
  maxoverlap_of_methods_by_sampling[[sampname]] <- sampmatrix_for_interM_interS
}

### INTER-METHOD INTER-SAMPLING
# Use overlaps max methods intra-sampling (cf ### INTER-METHOD INTRA-SAMPLING) to compare to other sampling
dir.create(file.path(outputpath, 'STATS/intermethod_intersampling'), showWarnings = FALSE)
interM_interS <- fromList(maxoverlap_of_methods_by_sampling)
interM_interS$intersect <- rowSums(interM_interS) 
write.csv(interM_interS,paste('STATS/intermethod_intersampling/MATRIX_intermethod_intersampling_feats_overlaps.csv',sep=""),row.names = TRUE)
upsetPlot <- upset(fromList(maxoverlap_of_methods_by_sampling),
                      sets = rev(names(maxoverlap_of_methods_by_sampling)),
                      keep.order = T,
                      nintersects = NA,
                      point.size = 9,
                      line.size = 2,
                      mb.ratio = c(0.48, 0.52),
                      matrix.dot.alpha = 0.6,
                      text.scale = 4)
ggsave(paste('STATS/intermethod_intersampling/UPSETPLOT_intermethod_intersampling_feats_overlaps.png',sep=""),ggplotify::as.ggplot(upsetPlot),width = 24,height = 9)
  
fileConn<-file("finished.txt")
writeLines(c("BDML to NEO4J and Statistics analysis is done!","Check STATS/ folder and intra-run folders."), fileConn)
close(fileConn)