#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Train meta-classifier to polish LOHGIC predictions for germline alleles in PROFILE


#########
# Setup #
#########
require(argparse, quietly=TRUE)
options(scipen=1000, stringsAsFactors=F)


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Train and apply classifier to polish",
                                           "LOHGIC germline predictions"))
parser$add_argument("--lohgic", metavar=".tsv", type="character", action="append",
                    help="Extended OncDRS mutations with LOHGIC predictions",
                    required=TRUE)
parser$add_argument("--invitae", metavar=".tsv", type="character",
                    help="Invitae panel data for both training & validation",
                    required=TRUE)
parser$add_argument("--seq-metrics", metavar=".tsv", type="character",
                    help="Panel sequencing metrics for each sample",
                    required=TRUE)
parser$add_argument("--known-somatic-drivers", metavar=".tsv", type="character",
                    help=paste("List of established cancer driver mutations. Two-column",
                               ".tsv of gene symbol and protein consequence (e.g., p.G12D)"),
                    required=TRUE)
parser$add_argument("--profile-frequencies", metavar=".tsv", type="character",
                    help=paste("Frequency per mutation in other patients in PROFILE"),
                    required=TRUE)
# parser$add_argument("--driver-universe", metavar=".txt", type="character",
#                     help="List of genes considered/tested for --known-somatic-drivers")
parser$add_argument("--validation-ids", metavar="string", type="character",
                    help=paste("Sample IDs to hold out from training for",
                               "validation; by default, all samples will be",
                               "used in training and validation will not be",
                               "conducted"))
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="output .tsv file for mutations with predictions")
args <- parser$parse_args()

# Load necessary libraries
require(OncoModR, quietly=TRUE)
require(caret, quietly=TRUE)
# require(randomForest, quietly=TRUE) # required for model = "rf"
# require(adabag, quietly=TRUE) # required for model = "AdaBag"
# require(fastAdaboost, quietly=TRUE) # required for model = "adaboost"
# require(RRF, quietly=TRUE) # required for model = "RRF"
# require(caTools, quietly=TRUE) # required for model = "LogitBoost"
require(ranger, quietly=TRUE) # required for model = "ranger"
# require(RWeka) # required for model = "LMT"
# require(stepPlr) #required for model = "plr"
require(plyr, quietly=TRUE)

# # DEV
# setwd("~/Desktop/Collins/VanAllen/PROFILE/LOHGIC_2023/")
# args <- list("lohgic" = "data/PROFILE.LOHGIC.annotated.tsv.gz",
#              "invitae" = "data/invitae_data.raw.txt.gz",
#              "seq_metrics" = "data/PROFILE.oncopanel_seq_metrics.subset.tsv.gz",
#              "known_somatic_drivers" = "data/msk_hotspots_v2.simple.tsv",
# #              "driver_universe" = "data/msk_impact_genes.list",
#              "profile_frequencies" = "data/PROFILE.mut_freqs.hotspot_test.tsv",
#              "validation_ids" = "data/LOHGIC.held_out_validation_samples.list",
#              "outfile" = "data/PROFILE.LOHGIC.annotated.final_predictions.tsv",
#              "feat_impt_pdf" = "plots/PROFILE.classifier.feature_importance.pdf")

# Read input data
muts <- read.table(args$lohgic, sep="\t", quote="", header=T, comment.char="")
seq.mx <- read.table(args$seq_metrics, sep="\t", quote="", header=T, comment.char="")
muts <- merge(muts, seq.mx, by="SAMPLE_ACCESSION_NBR", all.x=T, all.y=F, sort=F)
muts$PCT_TARGET_BASES_LT50X[which(is.na(muts$PCT_TARGET_BASES_LT50X))] <- 0
drivers <- read.table(args$known_somatic_drivers, header=T, sep="\t", comment.char="")
colnames(drivers) <- c("CANONICAL_GENE", "CANONICAL_PROTEIN_CHANGE")
drivers$known_somatic_driver <- TRUE
muts <- merge(muts, drivers, sort=F, all.x=T, all.y=F)
muts$known_somatic_driver[which(is.na(muts$known_somatic_driver))] <- FALSE
# driver.universe <- unique(read.table(args$driver_universe, header=F))
# colnames(driver.universe) <- "CANONICAL_GENE"
# driver.universe$on_impact_panel <- TRUE
# muts <- merge(muts, driver.universe, sort=F, all.x=T, all.y=F)
# muts$on_impact_panel[which(is.na(muts$on_impact_panel))] <- FALSE

idat <- read.table(args$invitae, sep="\t", quote="", header=T, comment.char="")
if(!is.null(args$validation_ids)){
  val.ids <- as.character(unique(read.table(args$validation_ids)[, 1]))
}else{
  val.ids <- c()
}

# Add derived features to mutation data as needed
muts$MAX_GNOMAD_FREQUENCY[which(is.na(muts$MAX_GNOMAD_FREQUENCY))] <- 0
muts$gnomad_overlap <- muts$MAX_GNOMAD_FREQUENCY > 0
muts$gnomad_common <- muts$MAX_GNOMAD_FREQUENCY > 0.001
muts$norm_depth <- muts$DEPTH / muts$MEAN_SAMPLE_COVERAGE
muts$DEL <- relevel(as.factor(muts$PLOIDY < 2), ref="FALSE")
muts$AMP <- relevel(as.factor(muts$PLOIDY > 2), ref="FALSE")
muts$model1 <- relevel(as.factor(sapply(muts$lohgic_model_1, function(s){
  unlist(strsplit(gsub(",", "", s), split=" CN_"))[1]
})), ref="Somatic")
muts$model2 <- relevel(as.factor(sapply(muts$lohgic_model_2, function(s){
  unlist(strsplit(gsub(",", "", s), split=" CN_"))[1]
})), ref="Somatic")
# muts$m1 <- sapply(muts$lohgic_model_1, function(s){unlist(strsplit(gsub(",", "", s), split=" "))[1]})
# muts$m2 <- sapply(muts$lohgic_model_2, function(s){unlist(strsplit(gsub(",", "", s), split=" "))[1]})
# muts$lohgic_concordance <- as.factor(paste(muts$m1, muts$m2, sep="_"))
muts$lohgic_mut_CN <- as.numeric(sapply(muts$lohgic_model_1, function(s){
  unlist(strsplit(s, split="{mut}=", fixed=T))[2]
}))
muts$gene_in_driver_set <- muts$CANONICAL_GENE %in% drivers$CANONICAL_GENE

# Annotate with frequency in other samples in PROFILE (matched on test type & panel version)
pro.freq <- read.table(args$profile_frequencies, header=T, sep="\t", comment.char="")
muts$CANONICAL_CHANGE <- muts$CANONICAL_PROTEIN_CHANGE
muts$CANONICAL_CHANGE[which(is.na(muts$CANONICAL_CHANGE))] <- muts$CANONICAL_CDNA_CHANGE[which(is.na(muts$CANONICAL_CHANGE))]
muts <- merge(muts, pro.freq, by=c("TEST_TYPE", "PANEL_VERSION", "CANONICAL_CHANGE", "CANONICAL_GENE"),
              all.x=T, all.y=F, sort=F)
muts$frac_other_patients_with_mut[which(is.na(muts$frac_other_patients_with_mut))] <- 0

# Define subset of samples to use for training
# Also restrict to only OncoPanel, not RapidHeme, run on solid tumor biopsies
exclude.biopsies <- c("BLOOD", "DIAGNOSIS_PENDING", "DIAGNOSIS_UNAVAILABL", "LYMPH",
                      "NEW_NODE_PENDING", "NO_ONCOTREE_NODE_FOU", "OTHER", "UNKNOWN")
all.ids <- unique(muts$SAMPLE_ACCESSION_NBR[which(muts$DFCI_MRN %in% idat$mrn & !(muts$BIOPSY_SITE %in% exclude.biopsies))])
cat(paste("Identified", prettyNum(length(all.ids), big.mark=","),
          "overlapping MRNs between --lohgic and --invitae inputs.\n"))
val.ids <- intersect(val.ids, all.ids)
# set.seed(2023)
# val.ids <- sample(all.ids, floor(0.1 * length(all.ids)), replace=F)
train.ids <- setdiff(all.ids, val.ids)
cat(paste("Dividing eligible samples into", prettyNum(length(train.ids), big.mark=","),
          "training samples and", prettyNum(length(val.ids), big.mark=","),
          "held-out validation samples.\n"))

# Merge data for training and validation
mdat <- merge(muts[which(muts$DFCI_MRN %in% idat$mrn), ],
              idat[which(idat$mrn %in% muts$DFCI_MRN), ],
              by.x=c("CHROMOSOME", "POSITION", "DFCI_MRN"),
              by.y=c("chrom", "vcf_pos", "mrn"),
              all.x=T, all.y=F, sort=F)

# Filter training.validation data to genes present in both OncoPanel and Invitae data
op.genes <- unique(muts$CANONICAL_GENE)
i.genes <- unique(idat$gene)
train.genes <- intersect(op.genes, i.genes)
cat(paste("Found", prettyNum(length(train.genes), big.mark=","),
          "gene symbols overlapping between --lohgic and --invitae inputs.\n"))
mdat <- mdat[which(mdat$CANONICAL_GENE %in% train.genes
                   & !(mdat$TEST_TYPE %in% "RAPIDHEME_CLINICAL")
                   & !(mdat$BIOPSY_SITE %in% exclude.biopsies)), ]

# Exclude genes from true negative set that have at least one high-confidence
# LOHGIC germline mutation but have <33% validation rate
# These are suspicious genes and shouldn't bias our results (e.g., TP53)
vrate.pergene <- as.data.frame(do.call("rbind", lapply(train.genes, function(g){
  idxs <- which(mdat$CANONICAL_GENE == g)
  pred.germ <- intersect(idxs, which(mdat$lohgic_germline_weight > 0.7))
  val.germ <- intersect(pred.germ, which(!is.na(mdat$hgvs)))
  c(g, length(idxs), length(pred.germ), length(val.germ), length(val.germ) / length(pred.germ))
})))
colnames(vrate.pergene) <- c("gene", "total_muts", "total_hq_germline", "hq_germline_val", "val_rate")
bad.genes <- vrate.pergene$gene[which(vrate.pergene$total_hq_germline > 0
                                      & vrate.pergene$val_rate < 1/3)]
cat(paste("Identified", prettyNum(length(bad.genes), big.mark=","),
          "suspicious genes with low germline validation rates:",
          paste(sort(bad.genes), collapse=", "), "\n"))

# Extract training and validation data from Invitae overlaps
features <- c("lohgic_germline_weight", "VAF", "norm_depth",
              "MEAN_SAMPLE_COVERAGE", "PURITY", "PURITY_CI",
              "gnomad_overlap", "MAX_GNOMAD_FREQUENCY",
              "known_somatic_driver", "gene_in_driver_set",
              "frac_other_patients_with_mut",
              "PLOIDY", "PANEL_VERSION", "model1", "model2",
              "lohgic_mut_CN", "SAMPLE_ACCESSION_NBR")
tp <- mdat[which(!is.na(mdat$hgvs)), features]
tp$VAL <- "germline"
tp$ground <- TRUE
tn <- mdat[which(is.na(mdat$hgvs) & !(mdat$gene %in% bad.genes)), features]
tn$VAL <- "not"
tn$ground <- TRUE

# For training, also add all variants that appear at AF>10% in gnomAD
# are not deemed pathogenic, and have VAF ~ [0.4, 0.6], as these are most likely germline
tp.snp.idx <- which(muts$MAX_GNOMAD_FREQUENCY >= 0.1
                    & muts$VAF >= 0.4
                    & muts$VAF <= 0.6
                    & muts$PATHOLOGIST_TIER > 3
                    & !(muts$SAMPLE_ACCESSION_NBR %in% c(train.ids, val.ids))
                    & muts$TEST_TYPE != "RAPIDHEME_CLINICAL"
                    & !(muts$BIOPSY_SITE %in% exclude.biopsies))
tp.snp <- muts[tp.snp.idx, features]
tp.snp$VAL <- "germline"
tp.snp$ground <- FALSE

# For training, also add 500 randomly selected activating tier 1 mutations
# that never appear in gnomAD and don't have VAF ~ [0.4, 0.6] as true negatives,
# as these are extremely unlikely to be germline
# Note that we restricted this specifically to activating mutations as het LoF of
# TSGs are much more common germline predisposition alleles, whereas that
# is not as often the case for oncogenes
tn.onc.idxs <- intersect(which(muts$PATHOLOGIST_TIER == 1
                               & !(muts$gnomad_overlap)
                               & (muts$VAF < 0.4 | muts$VAF > 0.6)
                               & !(muts$SAMPLE_ACCESSION_NBR %in% c(train.ids, val.ids))
                               & muts$TEST_TYPE != "RAPIDHEME_CLINICAL"
                               & !(muts$BIOPSY_SITE %in% exclude.biopsies)),
                         grep("missense", tolower(muts$CANONICAL_VARIANT_CLASS)))
set.seed(2023)
tn.onc <- muts[sample(tn.onc.idxs, 500), features]
tn.onc$VAL <- "not"
tn.onc$ground <- FALSE

# Combine + clean various subsets for training & testing
mdat.sub <- do.call("rbind", list(tp, tn, tp.snp, tn.onc))
mdat.sub$VAL <- relevel(as.factor(mdat.sub$VAL), ref="not")
mdat.sub$PANEL_VERSION <- relevel(as.factor(mdat.sub$PANEL_VERSION), ref="1")
tdat <- mdat.sub[which(!(mdat.sub$SAMPLE_ACCESSION_NBR %in% val.ids)), ]
tdat$SAMPLE_ACCESSION_NBR <- NULL
tdat$ground <- NULL
vdat <- mdat.sub[which(mdat.sub$SAMPLE_ACCESSION_NBR %in% val.ids & mdat.sub$ground), ]
vdat$ground <- NULL

# Train random forest using 10-fold cross-validation
cat(paste("Using a total of", prettyNum(length(which(tdat$VAL == "germline")), big.mark=","),
          "positive and", prettyNum(length(which(tdat$VAL == "not")), big.mark=","),
          "negative mutations for training.\n"))
trcontrol <- trainControl(method="cv", number=10, p=0.8,
                          classProbs = TRUE, summaryFunction = twoClassSummary,
                          savePredictions = "all")
set.seed(2023)
model <- train(VAL ~ . , data=tdat, method = "ranger",
               trControl = trcontrol, metric="ROC")
tdat.wpred <- cbind(tdat, predict(model, newdata=tdat, "prob"))
thresh.cut.steps <- seq(0, 1, 0.01)
thresh.cut <- as.data.frame(t(sapply(thresh.cut.steps, function(k){
  confusionMatrix(as.factor(sapply(tdat.wpred$germline, function(v){if(v > k){"germline"}else{"not"}})),
                  tdat.wpred$VAL, positive="germline")$byClass
})))
f.score <- 1
f.score.col <- paste("f", f.score, sep="")
thresh.cut[, f.score.col] <- (1+(f.score^2)) * ((thresh.cut$Precision*thresh.cut$Recall) / ((f.score^2 * thresh.cut$Precision) + thresh.cut$Recall))
threshold <- thresh.cut.steps[min(which(thresh.cut[, f.score.col] == max(thresh.cut[, f.score.col], na.rm=T)))]
cat(paste("Optimized decision threshold of P(Germline) > ", threshold,
          " using F_", f.score, " score on training data.\n", sep=""))
tdat.wpred$PRED <- sapply(tdat.wpred$germline, function(v){if(v > threshold){"germline"}else{"not"}})
tdat.wpred$PRED <- relevel(as.factor(tdat.wpred$PRED), ref="not")
cat("Finished training. Summary of training performance:\n")
confusionMatrix(tdat.wpred$PRED, tdat.wpred$VAL, positive="germline")

# Validate on held-out data, if optioned
if(nrow(vdat) > 0){
  vdat.wpred <- cbind(vdat, predict(model, newdata=vdat, "prob"))
  vdat.wpred$PRED <- sapply(vdat.wpred$germline, function(v){if(v>threshold){"germline"}else{"not"}})
  vdat.wpred$PRED <- relevel(as.factor(vdat.wpred$PRED), ref="not")
  cat("Summary of validation performance:\n")
  confusionMatrix(vdat.wpred$PRED, vdat.wpred$VAL, positive="germline")
}

# Compare results from using flat LOHGIC weight cutoffs of 0.5 and 0.7
for(cutoff in c(0.5, 0.7)){
  cat(paste("For comparison, summary of performance using fixed LOHGIC germline",
            "weight cutoff >", cutoff, ":\n"))
  print(confusionMatrix(relevel(as.factor(sapply(vdat.wpred$lohgic_germline_weight, function(v){if(v > cutoff){"germline"}else{"not"}})), ref="not"),
                        vdat.wpred$VAL, positive="germline"))
}

# Apply to all data
eligible.muts <- which(complete.cases(muts[, features]))
muts.complete <- muts[eligible.muts, ]
muts.complete$PANEL_VERSION <- relevel(as.factor(muts.complete$PANEL_VERSION), ref="1")
muts.complete$prob_germline <- predict(model, muts.complete, "prob")$germline
muts.complete$pred_label <- sapply(muts.complete$prob_germline, function(v){if(v >= threshold){"germline"}else{"somatic"}})
muts.all <- rbind.fill(muts.complete, muts[-eligible.muts, ])
write.table(muts.all, args$outfile, sep="\t", col.names=T, row.names=F, quote=F)

# # Plot feature importance, if optioned
# if(!is.null(args$feat_impt_pdf)){
#   pdf(args$`feat-impt-pdf`, height=4.5, width=5)
#   varImpPlot(model$finalModel, pch=19, main="Feature Importance", col="darkblue")
#   dev.off()
# }
