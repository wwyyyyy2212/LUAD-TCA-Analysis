rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(harmony)
  library(qs)
})

#------------------------------ Input/output -------------------------------
copykat_result_file <- "CopyKAT_prediction.txt"  # externally generated CopyKAT result
signature_file <- "TCA.txt"                 
out_dir <- "single_cell_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
scRNA<-readRDS('scRNA.rds')
scRNA <- subset(scRNA, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)

pdf(file.path(out_dir, "QC_after_filtering.pdf"), width = 9, height = 4)
print(VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

#------------------------------  Normalize/PCA -------------------------------
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA), verbose = FALSE)
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA), verbose = FALSE)

pdf(file.path(out_dir, "ElbowPlot.pdf"), width = 6, height = 4)
print(ElbowPlot(scRNA, ndims = 30))
dev.off()

#------------------------------  Harmony/clustering -------------------------
# Harmony is run only if a patient/batch column is present.
if ("patient" %in% colnames(scRNA@meta.data)) {
  scRNA <- RunHarmony(scRNA, group.by.vars = "patient", reduction = "pca", dims.use = 1:30)
  reduction_use <- "harmony"
} else {
  warning("No 'patient' metadata column found; using PCA rather than Harmony.")
  reduction_use <- "pca"
}

scRNA <- FindNeighbors(scRNA, reduction = reduction_use, dims = 1:15)
scRNA <- FindClusters(scRNA, resolution = 0.1)
scRNA <- RunUMAP(scRNA, reduction = reduction_use, dims = 1:15)

pdf(file.path(out_dir, "UMAP_clusters.pdf"), width = 7, height = 6)
print(DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
dev.off()

#------------------------------ Cluster markers ----------------------------
cluster_markers <- FindAllMarkers(
  scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
)
write.csv(cluster_markers, file.path(out_dir, "cluster_markers.csv"), row.names = FALSE)

#------------------------------ Major cell-type annotation -----------------
# Marker panel is based on the supplied script. Inspect DotPlot/markers and
# confirm cluster IDs before using the mapping below on a different dataset.
marker_panel <- list(
  T = c("CD3D", "CD3E"),
  B = c("CD79A", "MS4A1"),
  Plasma = c("TNFRSF17", "JCHAIN"),
  Myeloid = c("LYZ", "AIF1", "HLA-DQB2", "HLA-DPB1"),
  Epithelial = c("EPCAM", "KRT8"),
  Endothelial = c("PECAM1", "VWF"),
  Fibroblast = c("COL3A1", "DCN"),
  Mast = c("TPSAB1", "TPSB2")
)
marker_panel <- lapply(marker_panel, function(x) intersect(x, rownames(scRNA)))
pdf(file.path(out_dir, "Celltype_marker_DotPlot.pdf"), width = 10, height = 6)
print(DotPlot(scRNA, features = marker_panel) + RotatedAxis())
dev.off()
#------------------------------ Malignant-cell identification ---------------

if (!file.exists(copykat_result_file)) {
  warning("CopyKAT result file not found. Skipping malignant-cell annotation block.")
} else {
  copykat_pred <- read.table(copykat_result_file, header = TRUE, sep = "\t",
                             row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  pred_col <- if ("copykat.pred" %in% colnames(copykat_pred)) "copykat.pred" else NA_character_
  if (is.na(pred_col)) stop("CopyKAT file must contain a 'copykat.pred' column.")

  scRNA$copykat_status <- NA_character_
  common_cells <- intersect(colnames(scRNA), rownames(copykat_pred))
  scRNA$copykat_status[match(common_cells, colnames(scRNA))] <- as.character(copykat_pred[common_cells, pred_col])
  scRNA$copykat_status <- str_replace_all(scRNA$copykat_status, "not.defined", "diploid")
  scRNA$merge_celltype <- as.character(scRNA$celltype)
  scRNA$merge_celltype[scRNA$celltype == "Epithelial" & scRNA$copykat_status == "aneuploid"] <- "Malignant_Epithelial"
  scRNA$merge_celltype[scRNA$celltype == "Epithelial" & scRNA$copykat_status == "diploid"] <- "Non-malignant_Epithelial"

  write.csv(scRNA@meta.data, file.path(out_dir, "metadata_with_malignant_annotation.csv"))
  pdf(file.path(out_dir, "UMAP_malignant_annotation.pdf"), width = 7, height = 6)
  print(DimPlot(scRNA, reduction = "umap", group.by = "merge_celltype", label = FALSE))
  dev.off()
}

#------------------------------ Signature scoring ---------------------------
if (file.exists(signature_file)) {
  
  sig_df <- read.table(
    signature_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if ("genes" %in% colnames(sig_df)) {
    signature_genes <- unique(as.character(sig_df$genes))
  } else {
    stop("The signature file does not contain a column named 'genes'.")
  }
  
  signature_genes <- signature_genes[
    !is.na(signature_genes) & signature_genes != ""
  ]
  signature_genes <- intersect(signature_genes, rownames(scRNA))
  
  if (length(signature_genes) < 2) {
    
    warning("Fewer than two signature genes matched the expression matrix; scoring skipped.")
    
  } else {
    
    gene_sets <- list(gene = signature_genes)
    DefaultAssay(scRNA) <- "RNA"
    
    # AUCell scoring
    if (requireNamespace("AUCell", quietly = TRUE)) {
      
      expr_data <- GetAssayData(
        scRNA,
        assay = "RNA",
        slot = "data"
      )
      
      rankings <- AUCell::AUCell_buildRankings(
        expr_data,
        nCores = 1,
        plotStats = FALSE
      )
      
      auc_obj <- AUCell::AUCell_calcAUC(
        gene_sets,
        rankings,
        aucMaxRank = ceiling(nrow(rankings) * 0.1),
        nCores = 1
      )
      
      auc_mat <- AUCell::getAUC(auc_obj)
      
      scRNA$AUCell_gene <- as.numeric(
        auc_mat["gene", colnames(scRNA)]
      )
    }
    
    # UCell and singscore scoring
    if (requireNamespace("irGSEA", quietly = TRUE)) {
      
      scRNA <- irGSEA::irGSEA.score(
        object = scRNA,
        assay = "RNA",
        slot = "data",
        seeds = 123,
        ncores = 1,
        msigdb = FALSE,
        custom = TRUE,
        geneset = gene_sets,
        method = c("UCell", "singscore"),
        maxGSSize = 1000
      )
      
      if ("UCell" %in% names(scRNA@assays)) {
        
        ucell_mat <- GetAssayData(
          scRNA,
          assay = "UCell",
          slot = "counts"
        )
        
        if ("gene" %in% rownames(ucell_mat)) {
          scRNA$UCell_gene <- as.numeric(
            ucell_mat["gene", colnames(scRNA)]
          )
        }
      }
      
      if ("singscore" %in% names(scRNA@assays)) {
        
        sing_mat <- GetAssayData(
          scRNA,
          assay = "singscore",
          slot = "counts"
        )
        
        if ("gene" %in% rownames(sing_mat)) {
          scRNA$singscore_gene <- as.numeric(
            sing_mat["gene", colnames(scRNA)]
          )
        }
      }
    }
    
    # Seurat AddModuleScore
    scRNA <- AddModuleScore(
      scRNA,
      features = list(signature_genes),
      name = "AddModuleScore_gene"
    )
    
    if ("AddModuleScore_gene1" %in% colnames(scRNA@meta.data)) {
      scRNA$AddModuleScore_gene <- scRNA$AddModuleScore_gene1
    }
    
    # Integrate scores from different methods
    score_cols <- intersect(
      c(
        "AUCell_gene",
        "UCell_gene",
        "singscore_gene",
        "AddModuleScore_gene"
      ),
      colnames(scRNA@meta.data)
    )
    
    if (length(score_cols) >= 2) {
      
      score_df <- as.data.frame(
        scRNA@meta.data[, score_cols, drop = FALSE]
      )
      
      # Standardize scores
      score_scaled <- as.data.frame(scale(score_df))
      
      # Normalize standardized scores to the 0-1 range
      normalize01 <- function(x) {
        
        rng <- range(x, na.rm = TRUE)
        
        if (!all(is.finite(rng)) || diff(rng) == 0) {
          return(rep(0, length(x)))
        }
        
        (x - rng[1]) / diff(rng)
      }
      
      score_norm <- as.data.frame(
        lapply(score_scaled, normalize01)
      )
      
      colnames(score_norm) <- paste0(
        colnames(score_norm),
        "_scaled01"
      )
      
      # Calculate the composite score
      scRNA$Composite_Scoring <- rowSums(
        score_norm,
        na.rm = TRUE
      )
      
      # Export scoring results
      write.csv(
        cbind(
          cell = colnames(scRNA),
          scRNA@meta.data[
            ,
            c(score_cols, "Composite_Scoring"),
            drop = FALSE
          ]
        ),
        file = file.path(out_dir, "signature_scores.csv"),
        row.names = FALSE
      )
      
    } else {
      
      warning(
        "Fewer than two scoring methods produced results; composite scoring skipped."
      )
    }
  }
  
} else {
  
  warning("Signature gene file not found; scoring block skipped.")
  
}

#------------------------------ Malignant-cell high/low groups --------------
# The original workflow used the median of the composite score in malignant
# epithelial cells. This block runs only when both the annotation and score exist.
if ("merge_celltype" %in% colnames(scRNA@meta.data) && "Composite_Scoring" %in% colnames(scRNA@meta.data)) {
  malignant <- subset(scRNA, subset = merge_celltype == "Malignant_Epithelial")
  if (ncol(malignant) > 0) {
    score_cutoff <- median(malignant$Composite_Scoring, na.rm = TRUE)
    malignant$GT_group <- ifelse(malignant$Composite_Scoring > score_cutoff, "TCAhighepi", "CAlowepi")
    write.csv(malignant@meta.data, file.path(out_dir, "malignant_epithelial_high_low_metadata.csv"))

    if (length(unique(malignant$GT_group)) == 2) {
      Idents(malignant) <- "GT_group"
      deg <- FindMarkers(malignant, ident.1 = "TCAhighepi", ident.2 = "TCAlowepi",
                         logfc.threshold = 0.25, min.pct = 0.1)
      deg$gene <- rownames(deg)
      write.csv(deg, file.path(out_dir, "high_vs_low_DEGs.csv"), row.names = FALSE)
    }
    saveRDS(malignant, file.path(out_dir, "malignant_epithelial_high_low.rds"))
  }
}

#------------------------------ Save object -------------------------------
saveRDS(scRNA, file.path(out_dir, "seurat_processed_annotated.rds"))
qsave(scRNA, file.path(out_dir, "seurat_processed_annotated.qs"))
write.csv(scRNA@meta.data, file.path(out_dir, "final_cell_metadata.csv"))
####cox####
###############################################################################
# Cox proportional hazards prognostic model workflow
# Based on the supplied Cox/LASSO script:
# differential genes -> expression/survival merge -> univariable Cox
# -> LASSO-Cox with cross-validation -> multivariable Cox/stepwise selection
# -> risk score and median split -> KM -> time-dependent ROC
# -> external validation -> optional independent prognostic analysis.
# Check input formats, endpoint coding, time units, and cohort identifiers.
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(survminer)
  library(ggplot2)
  library(timeROC)
})

#------------------------------ Input and Output ------------------------------
deg_file <- file.path("single_cell_results", "high_vs_low_DEGs.csv")
# Expression matrix: rows = genes, columns = samples; raw TPM/normalized expression.
expression_file <- "data_tpm.txt"
# Clinical table: row names = sample IDs, must contain futime and fustat.
clinical_file <- "clinical_survival.txt"
# External validation table: sample rows, futime/fustat and expression columns for selected genes.
external_file <- "validation_GSE13213.txt"
out_dir <- "prognostic_model_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------ Prepare Expression and Survival Data ----------
deg <- read.csv(deg_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!("gene" %in% colnames(deg))) {
  stop("The single-cell DEG file must contain a column named 'gene'.")
}
candidate_genes <- unique(as.character(deg$gene))
candidate_genes <- candidate_genes[!is.na(candidate_genes) & candidate_genes != ""]

expr <- read.table(expression_file, sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE, quote = "", comment.char = "")
clin <- read.table(clinical_file, sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE, quote = "", comment.char = "")

if (!all(c("futime", "fustat") %in% colnames(clin))) {
  stop("Clinical table must contain 'futime' and 'fustat'. futime should be survival time; fustat: 1=event, 0=censored.")
}

expr <- log2(as.matrix(expr) + 1)
common_genes <- intersect(candidate_genes, rownames(expr))
if (length(common_genes) < 2) {
  stop("Fewer than two single-cell candidate genes match the bulk expression matrix. Check gene identifiers and capitalization.")
}
expr <- expr[common_genes, , drop = FALSE]
expr_sample <- as.data.frame(t(expr), check.names = FALSE)
expr_sample$sample_id <- rownames(expr_sample)
clin$sample_id <- rownames(clin)
train <- merge(clin, expr_sample, by = "sample_id")
rownames(train) <- train$sample_id
train$sample_id <- NULL

# Assumes input futime is in days and converts to years, as in the source script.
# If the clinical file is already in years, change this to train$futime <- as.numeric(train$futime).
train$futime <- as.numeric(train$futime) / 365
train$fustat <- as.numeric(as.character(train$fustat))
train <- train[is.finite(train$futime) & !is.na(train$fustat) & train$futime >= 0, , drop = FALSE]
train$futime[train$futime <= 0] <- 0.003

# Only retain genes available in the merged training set.
gene_candidates <- intersect(common_genes, colnames(train))
if (length(gene_candidates) < 2) stop("No usable candidate gene columns after merging clinical and expression data.")

#------------------------------ Univariable Cox Regression ------------------
uni_results <- lapply(gene_candidates, function(gene) {
  dat <- train[, c("futime", "fustat", gene), drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  if (nrow(dat) < 10 || length(unique(dat[[gene]])) < 2) return(NULL)
  fit <- try(coxph(as.formula(paste0("Surv(futime, fustat) ~ `", gene, "`")), data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  sm <- summary(fit)
  data.frame(gene = gene,
             HR = unname(sm$conf.int[1, "exp(coef)"]),
             HR.95L = unname(sm$conf.int[1, "lower .95"]),
             HR.95H = unname(sm$conf.int[1, "upper .95"]),
             pvalue = unname(sm$coefficients[1, "Pr(>|z|)"]),
             stringsAsFactors = FALSE)
})
uni_results <- do.call(rbind, uni_results)
if (is.null(uni_results) || nrow(uni_results) == 0) stop("No valid univariable Cox models were fitted.")
uni_results <- uni_results[order(uni_results$pvalue), ]
write.table(uni_results, file.path(out_dir, "univariable_Cox.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Filter genes using the source workflow's nominal p < 0.05 threshold.
selected_genes <- uni_results$gene[uni_results$pvalue < 0.05]
if (length(selected_genes) < 2) {
  warning("Fewer than two genes pass univariable p<0.05; using the top genes by p-value for LASSO as a fallback.")
  selected_genes <- head(uni_results$gene, min(20, nrow(uni_results)))
}

#------------------------------ LASSO-Cox with CV --------------------------
x <- as.matrix(train[, selected_genes, drop = FALSE])
mode(x) <- "numeric"
y <- Surv(train$futime, train$fustat)
keep <- apply(x, 2, function(z) all(is.finite(z)) && sd(z, na.rm = TRUE) > 0)
x <- x[, keep, drop = FALSE]
if (ncol(x) < 1) stop("No variable genes available for LASSO-Cox.")

lasso_fit <- glmnet(x = x, y = y, family = "cox", alpha = 1)
pdf(file.path(out_dir, "LASSO_coefficient_paths.pdf"), width = 7, height = 6)
plot(lasso_fit, xvar = "lambda", label = TRUE)
dev.off()

cv_fit <- cv.glmnet(x = x, y = y, family = "cox", alpha = 1, nfolds = min(10, nrow(x)))
pdf(file.path(out_dir, "LASSO_cross_validation.pdf"), width = 7, height = 6)
plot(cv_fit)
abline(v = log(c(cv_fit$lambda.min, cv_fit$lambda.1se)), lty = 2)
dev.off()

lasso_coef <- as.matrix(coef(cv_fit, s = "lambda.min"))
lasso_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
if (length(lasso_genes) == 0) stop("LASSO selected no genes at lambda.min; inspect CV curve or use lambda.1se/lambda.min as appropriate.")
write.csv(data.frame(gene = rownames(lasso_coef), coefficient = lasso_coef[, 1]),
          file.path(out_dir, "LASSO_coefficients.csv"), row.names = FALSE)

#------------------------------ Multivariable Cox --------------------------
# Retain the source workflow's stepwise multivariable Cox model after LASSO.
model_dat <- train[, c("futime", "fustat", lasso_genes), drop = FALSE]
model_dat <- model_dat[complete.cases(model_dat), , drop = FALSE]
cox_formula <- as.formula(paste("Surv(futime, fustat) ~", paste(sprintf("`%s`", lasso_genes), collapse = " + ")))
full_cox <- coxph(cox_formula, data = model_dat, x = TRUE)
final_cox <- try(step(full_cox, direction = "both", trace = 0), silent = TRUE)
if (inherits(final_cox, "try-error")) final_cox <- full_cox
final_summary <- summary(final_cox)
final_genes <- names(coef(final_cox))
coef_table <- data.frame(
  gene = final_genes,
  coefficient = unname(coef(final_cox)),
  HR = unname(final_summary$conf.int[, "exp(coef)"]),
  HR.95L = unname(final_summary$conf.int[, "lower .95"]),
  HR.95H = unname(final_summary$conf.int[, "upper .95"]),
  pvalue = unname(final_summary$coefficients[, "Pr(>|z|)"]),
  row.names = NULL
)
write.csv(coef_table, file.path(out_dir, "final_multivariable_Cox_model.csv"), row.names = FALSE)
saveRDS(final_cox, file.path(out_dir, "final_Cox_model.rds"))

#------------------------------ Risk score and KM --------------------------
# Cox type='risk' returns relative risk; median split is calculated within the
# training cohort and the training cutoff is then reused for external validation.
train_score <- as.numeric(predict(final_cox, newdata = train, type = "risk"))
train_cutoff <- median(train_score, na.rm = TRUE)
train$riskScore <- train_score
train$risk <- ifelse(train$riskScore > train_cutoff, "high", "low")
train$risk <- factor(train$risk, levels = c("low", "high"))
write.table(cbind(sample_id = rownames(train), train[, c("futime", "fustat", final_genes, "riskScore", "risk"), drop = FALSE]),
            file.path(out_dir, "risk_training.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

km_fit <- survfit(Surv(futime, fustat) ~ risk, data = train)
km_plot <- ggsurvplot(km_fit, data = train, conf.int = TRUE, pval = TRUE,
                      risk.table = TRUE, risk.table.height = 0.25,
                      legend.title = "Risk group", legend.labs = c("Low risk", "High risk"),
                      xlab = "Time (years)", ylab = "Overall survival probability")
pdf(file.path(out_dir, "KM_training.pdf"), width = 7.5, height = 6)
print(km_plot)
dev.off()

# Risk score association with survival (univariable Cox).
risk_cox <- coxph(Surv(futime, fustat) ~ riskScore, data = train)
capture.output(summary(risk_cox), file = file.path(out_dir, "riskScore_Cox_summary.txt"))

#------------------------------ Time-dependent ROC -------------------------
roc_train <- timeROC(T = train$futime, delta = train$fustat, marker = train$riskScore,
                     cause = 1, weighting = "marginal", times = c(1, 3, 5), ROC = TRUE)
write.csv(data.frame(time_years = c(1, 3, 5), AUC = as.numeric(roc_train$AUC)),
          file.path(out_dir, "time_dependent_AUC_training.csv"), row.names = FALSE)
pdf(file.path(out_dir, "time_dependent_ROC_training.pdf"), width = 5, height = 5)
plot(roc_train, time = 1, title = FALSE)
plot(roc_train, time = 3, add = TRUE, title = FALSE)
plot(roc_train, time = 5, add = TRUE, title = FALSE)
legend("bottomright", legend = paste0(c(1, 3, 5), "-year AUC: ", sprintf("%.3f", roc_train$AUC)), bty = "n")
dev.off()

#------------------------------ External validation ------------------------
if (file.exists(external_file)) {
  ext <- read.table(external_file, sep = "\t", header = TRUE, row.names = 1,
                    check.names = FALSE, quote = "", comment.char = "")
  if (all(c("futime", "fustat", final_genes) %in% colnames(ext))) {
    ext$futime <- as.numeric(ext$futime) / 365  # change if already in years
    ext$fustat <- as.numeric(as.character(ext$fustat))
    ext$futime[ext$futime <= 0] <- 0.003
    ext <- ext[complete.cases(ext[, c("futime", "fustat", final_genes), drop = FALSE]), , drop = FALSE]
    ext_score <- as.numeric(predict(final_cox, newdata = ext, type = "risk"))
    ext$riskScore <- ext_score
    ext$risk <- factor(ifelse(ext_score > train_cutoff, "high", "low"), levels = c("low", "high"))
    write.table(cbind(sample_id = rownames(ext), ext[, c("futime", "fustat", final_genes, "riskScore", "risk"), drop = FALSE]),
                file.path(out_dir, "risk_external_validation.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    ext_km <- survfit(Surv(futime, fustat) ~ risk, data = ext)
    ext_plot <- ggsurvplot(ext_km, data = ext, conf.int = TRUE, pval = TRUE,
                           risk.table = TRUE, risk.table.height = 0.25,
                           legend.title = "Risk group", legend.labs = c("Low risk", "High risk"),
                           xlab = "Time (years)", ylab = "Overall survival probability")
    pdf(file.path(out_dir, "KM_external_validation.pdf"), width = 7.5, height = 6)
    print(ext_plot)
    dev.off()
    
    if (length(unique(ext$fustat)) > 1 && max(ext$futime, na.rm = TRUE) >= 5) {
      roc_ext <- timeROC(T = ext$futime, delta = ext$fustat, marker = ext$riskScore,
                         cause = 1, weighting = "marginal", times = c(1, 3, 5), ROC = TRUE)
      write.csv(data.frame(time_years = c(1, 3, 5), AUC = as.numeric(roc_ext$AUC)),
                file.path(out_dir, "time_dependent_AUC_external.csv"), row.names = FALSE)
      pdf(file.path(out_dir, "time_dependent_ROC_external.pdf"), width = 5, height = 5)
      plot(roc_ext, time = 1, title = FALSE)
      plot(roc_ext, time = 3, add = TRUE, title = FALSE)
      plot(roc_ext, time = 5, add = TRUE, title = FALSE)
      legend("bottomright", legend = paste0(c(1, 3, 5), "-year AUC: ", sprintf("%.3f", roc_ext$AUC)), bty = "n")
      dev.off()
    }
  } else {
    warning("External validation file lacks one or more required columns; external validation skipped.")
  }
}

#------------------------------ Independent prognostic value (optional) -----
clinical_vars <- setdiff(colnames(clin), c("futime", "fustat", "sample_id"))
if (length(clinical_vars) > 0) {
  clinical_for_model <- clin
  clinical_for_model$sample_id <- rownames(clin)
  common_ids <- intersect(rownames(train), clinical_for_model$sample_id)
  clin2 <- clinical_for_model[match(common_ids, clinical_for_model$sample_id), , drop = FALSE]
  risk2 <- train[common_ids, c("futime", "fustat", "riskScore"), drop = FALSE]
  independent_dat <- cbind(risk2, clin2[, clinical_vars, drop = FALSE])
  independent_dat <- independent_dat[complete.cases(independent_dat), , drop = FALSE]
  if (nrow(independent_dat) > 10) {
    uni_cov <- lapply(c("riskScore", clinical_vars), function(v) {
      if (length(unique(independent_dat[[v]])) < 2) return(NULL)
      fit <- try(coxph(as.formula(paste0("Surv(futime, fustat) ~ `", v, "`")), data = independent_dat), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      sm <- summary(fit)
      data.frame(variable = v, HR = sm$conf.int[1, "exp(coef)"],
                 HR.95L = sm$conf.int[1, "lower .95"], HR.95H = sm$conf.int[1, "upper .95"],
                 pvalue = sm$coefficients[1, "Pr(>|z|)"], row.names = NULL)
    })
    uni_cov <- do.call(rbind, uni_cov)
    if (!is.null(uni_cov)) write.csv(uni_cov, file.path(out_dir, "clinical_univariable_Cox.csv"), row.names = FALSE)
    
    # Multivariable model with riskScore and available clinical covariates.
    multi_vars <- c("riskScore", clinical_vars)
    multi_formula <- as.formula(paste("Surv(futime, fustat) ~", paste(sprintf("`%s`", multi_vars), collapse = " + ")))
    multi_fit <- try(coxph(multi_formula, data = independent_dat), silent = TRUE)
    if (!inherits(multi_fit, "try-error")) capture.output(summary(multi_fit), file = file.path(out_dir, "clinical_multivariable_Cox_summary.txt"))
  }
}

message("Prognostic model workflow completed. Confirm survival endpoint coding, time units, gene IDs, and external cohort expression scale before interpreting results.")
