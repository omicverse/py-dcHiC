#!/usr/bin/env Rscript
# dcHiC reference runner — invoked under conda env $R_TEST_ENV (R 4.x + dcHiC deps).
#
# Usage: Rscript r_reference_driver.R <fixture.json> <output.json>
#
# Computes the two scoped dcHiC outputs using the REAL upstream primitives:
#   Stage A (compartment calling): functionsdchic::oe2cor (C++) + svd + GC orientation
#   Stage B (differential):        limma::normalizeQuantiles + verbatim calcen +
#                                   robust::covRob + stats::mahalanobis + BH
# JSON keys match manifest.yaml::outputs[*].location_reference.

suppressMessages({
  library(jsonlite)
  library(functionsdchic)
  library(limma)
  library(robust)
})

args <- commandArgs(trailingOnly = TRUE)
fixture_path <- args[1]
output_path  <- args[2]
if (is.na(fixture_path) || is.na(output_path)) {
  stop("Usage: Rscript r_reference_driver.R <fixture.json> <output.json>")
}
fix <- fromJSON(fixture_path)

## ============================ Stage A ====================================== ##
A   <- fix$stageA
n   <- as.integer(A$n_bins)
res <- as.numeric(A$resolution)
a_idx  <- as.integer(A$a_idx)
b_idx  <- as.integer(A$b_idx)
weight <- as.numeric(A$weight)
pos    <- as.numeric(A$pos)
gcc    <- as.numeric(A$gcc)
count_thr <- 0; minexpcc <- 0; n_pcs <- 2

# Observed/Expected (dchicf.r:357-375)
kw <- weight > count_thr
a_idx <- a_idx[kw]; b_idx <- b_idx[kw]; weight <- weight[kw]
dst <- abs(pos[a_idx] - pos[b_idx])
agg <- aggregate(weight ~ dst, FUN = sum)               # cols: dst, weight
total <- sapply(agg$dst, function(d) { v <- n - d / res; sum(v[v > 0]) })
expcc <- agg$weight / total
floorv <- ifelse(min(expcc) > minexpcc, min(expcc), minexpcc)
expcc[expcc <= floorv] <- floorv
expcc_by_dist <- setNames(expcc, as.character(agg$dst))
weight_oe <- weight / expcc_by_dist[as.character(dst)]

mat  <- functionsdchic::ijk2mat(cbind(a_idx, b_idx, as.numeric(weight_oe)), n, n)
keep <- which(rowSums(mat) >= 3)
m    <- mat[keep, keep]
bins <- nrow(m)

# C1 = column correlation of O/E (real C++ oe2cor); mat2fbm finalize (diag=1, round5)
x   <- functionsdchic::oe2cor(m, c(0), c(bins - 1), 1, 0)
C1  <- x[[1]]$mat
diag(C1) <- 1; C1 <- round(C1, 5)
# C2 = column correlation of C1; keep z-scored C1 (zm2) for the projection
x2  <- functionsdchic::oe2cor(C1, c(0), c(bins - 1), 1, 0)
zm2 <- x2[[length(x2)]]$zmat
C2  <- x2[[1]]$mat
diag(C2) <- 1; C2 <- round(C2, 5)

sv <- svd(C2)
V  <- sv$v[, 1:n_pcs]
pc <- zm2 %*% V                                          # dchicf.r:300
gcc_keep <- gcc[keep]
for (k in 1:n_pcs) pc[, k] <- sign(cor(pc[, k], gcc_keep)) * pc[, k]
compartment_pc1 <- as.numeric(pc[, 1])

## ============================ Stage B ====================================== ##
# verbatim calcen (dchicf.r:751)
calcen <- function(df, class, rzscore, szscore) {
  df_dist <- as.data.frame(t(apply(df, 1, function(x) {
    sqrt(colSums(as.matrix(dist(as.numeric(x)))) / (length(as.numeric(x)) - 1))
  })))
  if (class == "rep") {
    value <- as.numeric(as.matrix(df_dist)); col_mean <- mean(value); col_sd <- sd(value)
    df_zsc <- list()
    for (nn in 1:ncol(df_dist)) df_zsc[[nn]] <- (df_dist[, nn] - col_mean) / col_sd
    df_zsc <- do.call(cbind, df_zsc)
    df_pvl <- round(pnorm(df_zsc, mean = rzscore, lower.tail = T), 5)
    df_pvl_max <- apply(df_pvl, 1, max)
    return(round(df * (1 - df_pvl_max), 5))
  } else {
    df_zsc <- list()
    for (nn in 1:ncol(df)) {
      value <- as.numeric(as.matrix(df_dist[, nn]))
      df_zsc[[nn]] <- (value - mean(value)) / sd(value)
    }
    df_zsc <- do.call(cbind, df_zsc)
    df_pvl <- round(pnorm(as.matrix(df_zsc), mean = szscore, lower.tail = T), 5)
    df_pvl_max <- apply(df_pvl, 1, max)
    return(round(df * (1 - df_pvl_max), 5))
  }
}

B <- fix$stageB
pc_raw <- as.matrix(B$pc_raw)                            # bins x replicates
storage.mode(pc_raw) <- "numeric"
conditions <- as.character(B$conditions)
uniq <- unique(conditions)
rzscore <- 2; szscore <- 0; refine <- TRUE; rconf <- 0.90

# per-condition quantile normalisation + mean-abs scaling (dchicf.r:880-882)
pc_qnm <- matrix(0, nrow(pc_raw), ncol(pc_raw))
for (c in uniq) {
  cols <- which(conditions == c)
  qn <- limma::normalizeQuantiles(pc_raw[, cols], ties = T)
  qn <- qn / apply(abs(qn), 2, mean)
  pc_qnm[, cols] <- qn
}

cen_list <- list(); grp_list <- list()
for (c in uniq) {
  cols <- which(conditions == c)
  sub  <- pc_qnm[, cols, drop = FALSE]
  if (ncol(sub) > 1) {
    cen_list[[c]] <- as.matrix(calcen(as.data.frame(sub), "rep", rzscore, szscore))
    grp_list[[c]] <- apply(sub, 1, mean)
  } else {
    cen_list[[c]] <- rep(mean(sub), nrow(sub))
    grp_list[[c]] <- sub[, 1]
  }
}
intra_grp     <- do.call(cbind, grp_list)
ncond         <- ncol(intra_grp)
intra_grp_cen <- as.matrix(calcen(as.data.frame(intra_grp), "sample", rzscore, szscore))
max_min_diff  <- abs(apply(intra_grp, 1, max) - apply(intra_grp, 1, min))^2

maha_pval <- function(cov_rows) {
  sam_cov <- robust::covRob(cov_rows)$cov
  inv_cov <- solve(sam_cov)
  smaha <- numeric(nrow(intra_grp))
  for (k in 1:nrow(intra_grp)) {
    smaha[k] <- mahalanobis(intra_grp[k, ], as.matrix(intra_grp_cen[k, ]),
                            (inv_cov * max_min_diff[k]), inverted = TRUE)
  }
  list(smaha = smaha, pval = pchisq(smaha, df = ncond - 1, lower.tail = FALSE))
}

r   <- maha_pval(intra_grp)
pv  <- r$pval
if (refine) {
  thr <- pchisq(qchisq(rconf, df = ncond - 1), df = ncond - 1, lower.tail = FALSE)
  r   <- maha_pval(intra_grp[pv > thr, ])
  pv  <- r$pval
}
padj <- p.adjust(pv, "BH")

## ============================ write ======================================== ##
out <- list(
  compartment_pc1   = compartment_pc1,
  differential_padj = as.numeric(padj)
)
write_json(out, output_path, auto_unbox = TRUE, digits = NA, na = "null")
cat("[ref] wrote:", output_path, "\n")
