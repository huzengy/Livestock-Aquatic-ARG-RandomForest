# Random Forest validation: repeated cross-validation, importance agreement, and collinearity
# Run from the directory containing ARG_subtype_abundance.xlsx.
# Dependencies: readxl, randomForest
suppressPackageStartupMessages({library(readxl); library(randomForest)})
file_path <- "ARG_subtype_abundance.xlsx"
out_dir <- "results/validation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

water <- read_excel(file_path, sheet = "Water")
refs <- read_excel(file_path, sheet = "Sheet1")
label_col <- names(refs)[1]
water_X <- as.data.frame(water[-c(1, 2)])
river_X <- as.data.frame(refs[refs[[label_col]] == "Yangtze River", -c(1, 2)])
# Align common numeric ARG subtype columns across the two datasets.
common <- intersect(names(water_X)[vapply(water_X, is.numeric, logical(1))],
                    names(river_X)[vapply(river_X, is.numeric, logical(1))])
water_X <- water_X[, common, drop=FALSE]
river_X <- river_X[, common, drop=FALSE]
water_X[is.na(water_X)] <- 0
river_X[is.na(river_X)] <- 0
X <- rbind(water_X, river_X)
y <- factor(c(rep("Chongming", nrow(water_X)), rep("Yangtze", nrow(river_X))),
            levels=c("Chongming", "Yangtze"))
cat("Input samples: Chongming =", sum(y=="Chongming"), "; Yangtze =", sum(y=="Yangtze"), "\n")
if (sum(y=="Chongming") != 30 || sum(y=="Yangtze") != 23) {
  stop("Sample counts differ from manuscript (30/23). Verify input file before proceeding.")
}

# Apply the same directional enrichment criteria as the feature-selection script.
# During cross-validation, selection is repeated using training data only.
select_candidates <- function(x, labels) {
  a <- x[labels=="Chongming",,drop=FALSE]
  b <- x[labels=="Yangtze",,drop=FALSE]
  nms <- names(x)
  p <- vapply(nms, function(f) {
    tryCatch(suppressWarnings(wilcox.test(a[[f]], b[[f]], exact=FALSE)$p.value),
             error=function(e) 1)
  }, numeric(1))
  ma <- colMeans(a); mb <- colMeans(b)
  pa <- colMeans(a>0); pb <- colMeans(b>0)
  keep <- ma>mb & (ma+1e-10)/(mb+1e-10)>=1.5 & pa>pb & p.adjust(p,"BH")<0.05
  nms[keep]
}

score_metrics <- function(truth, prob, threshold=0.5) {
  predicted <- ifelse(prob >= threshold, "Chongming", "Yangtze")
  tab <- table(factor(truth,levels=levels(y)), factor(predicted,levels=levels(y)))
  sens <- if(sum(tab[1,])>0) tab[1,1]/sum(tab[1,]) else NA_real_
  spec <- if(sum(tab[2,])>0) tab[2,2]/sum(tab[2,]) else NA_real_
  # rank-based AUC (Chongming positive)
  positive <- truth == "Chongming"
  ranks <- rank(prob, ties.method="average")
  n1 <- sum(positive); n0 <- sum(!positive)
  auc <- (sum(ranks[positive]) - n1*(n1+1)/2)/(n1*n0)
  c(Accuracy=mean(predicted==truth), BalancedAccuracy=mean(c(sens,spec)), AUC=auc,
    Sensitivity=sens, Specificity=spec)
}

# Fit the full-data model using the same seed and predictor order as the primary analysis.
candidates <- select_candidates(X,y)
cat("Full-data directionally enriched candidate ARGs:",length(candidates),"\n")
if (length(candidates)<2) stop("Too few candidate ARGs; check input data compatibility")
Xcand <- X[,candidates,drop=FALSE]
Xcand <- Xcand[,colSums(Xcand)>0,drop=FALSE]
set.seed(123)
rf <- randomForest(x=Xcand,y=y,ntree=1000,importance=TRUE)
gini <- importance(rf,type=2)[,1]
perm <- importance(rf,type=1,scale=FALSE)[,1]
imp <- data.frame(ARG=names(gini),Gini=as.numeric(gini),Permutation=as.numeric(perm),stringsAsFactors=FALSE)
imp$Rank_Gini <- rank(-imp$Gini,ties.method="min")
imp$Rank_Permutation <- rank(-imp$Permutation,ties.method="min")
imp <- imp[order(imp$Rank_Gini),]
imp$Top20_Gini <- imp$Rank_Gini<=20
imp$Top20_Permutation <- imp$Rank_Permutation<=20
write.csv(imp,file.path(out_dir,"full_RF_feature_importance.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(as.data.frame.matrix(rf$confusion),file.path(out_dir,"OOB_confusion_matrix.csv"),fileEncoding="UTF-8")
oob <- data.frame(Method="OOB (after full-data filtering)",Error=tail(rf$err.rate[,"OOB"],1))
write.csv(oob,file.path(out_dir,"OOB_error.csv"),row.names=FALSE)

# Collinearity: rank correlations among full-data candidate features.
rho <- suppressWarnings(cor(Xcand,method="spearman",use="pairwise.complete.obs"))
diag(rho) <- NA_real_
high_pairs <- which(upper.tri(rho) & abs(rho)>=0.8,arr.ind=TRUE)
pairs <- if(nrow(high_pairs)) data.frame(ARG1=colnames(rho)[high_pairs[,1]],
  ARG2=colnames(rho)[high_pairs[,2]],rho=rho[high_pairs]) else data.frame(ARG1=character(),ARG2=character(),rho=numeric())
write.csv(pairs,file.path(out_dir,"correlated_ARG_pairs_absrho_ge_0.8.csv"),row.names=FALSE,fileEncoding="UTF-8")

# Repeated stratified 5-fold CV, all feature filtering done INSIDE each training fold.
# Five repetitions of stratified five-fold cross-validation (25 training/test splits).
REPEATS <- 5L
FOLDS <- 5L
set.seed(20261011)
cv_rows <- list(); cv_predictions <- list(); top20_freq <- setNames(integer(length(names(X))), names(X))
for (rep in seq_len(REPEATS)) {
  fold <- integer(length(y))
  for (cl in levels(y)) fold[which(y==cl)] <- sample(rep(seq_len(FOLDS),length.out=sum(y==cl)))
  p_all <- rep(NA_real_,length(y))
  for (k in seq_len(FOLDS)) {
    tr <- which(fold!=k); te <- which(fold==k)
    cc <- select_candidates(X[tr,,drop=FALSE],y[tr])
    if (length(cc)==0) stop(sprintf("No candidate passed filtering at CV repeat %d fold %d. Do NOT relax cutoffs silently.",rep,k))
    set.seed(2000 + 100*rep + k)
    m <- randomForest(x=X[tr,cc,drop=FALSE],y=y[tr],ntree=1000,importance=TRUE)
    vi <- sort(importance(m,type=2)[,1],decreasing=TRUE)
    selected <- names(head(vi,20))
    top20_freq[selected] <- top20_freq[selected]+1L
    # Refit on training fold's Top 20 to mirror final feature-selection target.
    set.seed(3000 + 100*rep + k)
    m20 <- randomForest(x=X[tr,selected,drop=FALSE],y=y[tr],ntree=1000)
    p_all[te] <- predict(m20,newdata=X[te,selected,drop=FALSE],type="prob")[,"Chongming"]
  }
  if(anyNA(p_all)) stop("CV prediction missing")
  cv_rows[[rep]] <- c(Repeat=rep,score_metrics(y,p_all))
  cv_predictions[[rep]] <- data.frame(Repeat=rep,Sample=seq_along(y),Class=y,P_Chongming=p_all)
  cat("Completed CV repetition",rep,"of",REPEATS,"\n")
}
cv_summary <- as.data.frame(do.call(rbind,cv_rows))
write.csv(cv_summary,file.path(out_dir,"cross_validation_metrics_by_repeat.csv"),row.names=FALSE)
write.csv(do.call(rbind,cv_predictions),file.path(out_dir,"cross_validation_out_of_fold_predictions.csv"),row.names=FALSE)
mean_sd <- function(v) sprintf("%.3f +/- %.3f",mean(v),sd(v))
cat("\nCV summary (mean +/- SD across repeats):\n")
for (col in c("Accuracy","BalancedAccuracy","AUC")) cat(col,mean_sd(cv_summary[[col]]),"\n")
topstab <- data.frame(ARG=names(top20_freq),Selected_Folds=as.integer(top20_freq),
                     Selection_Frequency=as.integer(top20_freq)/(REPEATS*FOLDS))
topstab <- topstab[order(-topstab$Selection_Frequency),]
write.csv(topstab,file.path(out_dir,"CV_top20_selection_frequency.csv"),row.names=FALSE,fileEncoding="UTF-8")

overlap <- sum(imp$Top20_Gini & imp$Top20_Permutation)
rank_rho <- suppressWarnings(cor(imp$Gini,imp$Permutation,method="spearman"))
cat("Gini/permutation Top 20 overlap:",overlap,"/20; Spearman rho of scores:",round(rank_rho,3),"\n")
cat("High-correlation candidate pairs (|rho| >= 0.8):",nrow(pairs),"\n")

cat("\nValidation results saved to:", normalizePath(out_dir), "\n")
