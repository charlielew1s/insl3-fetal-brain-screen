suppressPackageStartupMessages({
  library(ExpressionAtlas)
  library(SummarizedExperiment)
  library(edgeR)
  library(GenomeInfoDb)
})

dir.create("results", showWarnings=FALSE)
dir.create("logs", showWarnings=FALSE)

acc <- "E-MTAB-4840"
cat("Retrieving", acc, "from EMBL-EBI Expression Atlas...\n")
obj <- getAtlasExperiment(acc)

# Resolve RNA-seq SummarizedExperiment robustly.
if (methods::is(obj, "SummarizedExperiment")) {
  se <- obj
} else if (is.list(obj)) || methods::is(obj, "List")) {
  if ("rnaseq" %in% names(obj)) se <- obj[["rnaseq"]]
  else {
    ii <- which(vapply(obj, function(x) methods::is(x, "SummarizedExperiment"), logical(1)))
    if (!length(ii)) stop("No SummarizedExperiment found in getAtlasExperiment() result.")
    se <- obj[[ii[1]]]
  }
} else stop("Unexpected ExpressionAtlas object type: ", class(obj)[1])

assay_names <- SummarizedExperiment::assayNames(se)
cat("Assays:", paste(assay_names, collapse=", "), "\n")
count_name <- if ("counts" %in% assay_names) "counts" else assay_names[grep("count", assay_names, ignore.case=TRUE)[1]]
if (is.na(count_name) || is.null(count_name)) stop("No count assay found.")
cts <- SummarizedExperiment::assay(se, count_name)

meta <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors=FALSE)
meta$sample_id <- rownames(meta)

g <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors=FALSE)
g$gene_id <- rownames(g)

chr <- rep(NA_character_, nrow(g))
try({
  rr <- SummarizedExperiment::rowRanges(se)
  if (length(rr) == nrow(g)) chr <- as.character(GenomeInfoDb::seqnames(rr))
}, silent=TRUE)
if (all(is.na(chr))) {
  cc <- names(g)[grepl("chrom|seqname", names(g), ignore.case=TRUE)]
  if (length(cc)) chr <- as.character(g[[cc[1]]])
}
g$chromosome <- sub("^chr", "", chr, ignore.case=TRUE)

if (all(is.na(g$chromosome))) {
  stop("Could not resolve chromosome annotation. Refusing to fake an autosomal screen.")
}

# Normalize counts.
dge <- edgeR::DGEList(counts=cts)
keep <- edgeR::filterByExpr(dge)
dge <- dge[keep,, keep.lib.sizes=FALSE]
dge <- edgeR::calcNormFactors(dge, method="TMM")
lcpm <- edgeR::cpm(dge, log=TRUE, prior.count=0.5)
g <- g[match(rownames(lcpm), g$gene_id), , drop=FALSE]

# Metadata column discovery.
pick_col <- function(patterns, nms) {
  for (p in patterns) {
    h <- grep(p, nms, ignore.case=TRUE, value=TRUE)
    if (length(h)) return(h[1])
  }
  NA_character_
}
sex_col <- pick_col(c("^sex$", "karyotype", "fetal.*sex", "sex"), names(meta))
age_col <- pick_col(c("developmental.*stage", "post.*conception", "gestational.*age", "^age$", "age"), names(meta))
donor_col <- pick_col(c("individual", "donor", "fetus", "subject"), names(meta))
tissue_col <- pick_col(c("organism.*part", "brain.*region", "region", "tissue"), names(meta))

cat("Resolved columns:\n",
    "sex=", sex_col, "\nage=", age_col, "\ndonor=", donor_col, "\ntissue=", tissue_col, "\n", sep="")
if (any(is.na(c(sex_col,age_col,donor_col,tissue_col)))) {
  write.table(meta, "logs/metadata_debug.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  stop("Could not resolve essential metadata columns. See logs/metadata_debug.tsv")
}

norm_sex <- function(x) {
  z <- toupper(gsub("\\s+", "", as.character(x)))
  out <- rep(NA_character_, length(z))
  out[z %in% c("M","MALE","XY","46,XY","46XY")] <- "XY"
  out[z %in% c("F","FEMALE","XX","46,XX","46XX")] <- "XX"
  out[grepl("46,XY", z, fixed=TRUE)] <- "XY"
  out[grepl("46,XX", z, fixed=TRUE)] <- "XX"
  out[is.na(out) & grepl("MALE", z) & !grepl("FEMALE", z)] <- "XY"
  out[is.na(out) & grepl("FEMALE", z)] <- "XX"
  out
}

parse_age <- function(x) {
  z <- tolower(as.character(x))
  out <- rep(NA_real_, length(z))
  cs <- grepl("carnegie|\\bcs\\s*22|\\bcs\\s*23", z)
  out[cs] <- 7.75
  for (i in which(!cs)) {
    nums <- regmatches(z[i], gregexpr("[0-9]+(?:\\.[0-9]+)?", z[i], perl=TRUE))[[1]]
    nums <- suppressWarnings(as.numeric(nums))
    nums <- nums[is.finite(nums)]
    if (!length(nums)) next
    if (grepl("-", z[i], fixed=TRUE) && length(nums) >= 2) out[i] <- mean(nums[1:2])
    else out[i] <- nums[1]
  }
  out
}

meta$sex_std <- norm_sex(meta[[sex_col]])
meta$age_wpc <- parse_age(meta[[age_col]])
meta$donor_id <- as.character(meta[[donor_col]])
meta$region <- as.character(meta[[tissue_col]])

write.table(meta, "results/stage1_metadata_resolved.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Align sample IDs.
common <- intersect(meta$sample_id, colnames(lcpm))
if (!length(common)) stop("No overlap between colData rownames and expression columns.")
meta <- meta[match(common, meta$sample_id), , drop=FALSE]
lcpm <- lcpm[, common, drop=FALSE]

# Autosomal gene filtering.
chr_by_gene <- setNames(toupper(gsub("^CHR","", g$chromosome)), g$gene_id)
genes <- rownames(lcpm)
autosomal <- genes[!chr_by_gene[genes] %in% c("X","Y") & !is.na(chr_by_gene[genes])]
lcpm <- lcpm[autosomal,,drop=FALSE]
cat("Autosomal expressed genes:", nrow(lcpm), "\n")

# Keep resolved XX/XY fetal brain samples in 7-20 wpc.
use <- which(meta$sex_std %in% c("XX","XY") & is.finite(meta$age_wpc) &
             meta$age_wpc >= 7 & meta$age_wpc <= 20 &
             !is.na(meta$donor_id) & !is.na(meta$region))
meta <- meta[use,,drop=FALSE]
lcpm <- lcpm[,use,drop=FALSE]

# Collapse highly granular region labels to a reproducible cleaned label.
meta$region_clean <- trimws(tolower(meta$region))
meta$region_clean <- gsub("\\s+", " ", meta$region_clean)
meta$is_cortex_broad <- grepl("cortex|cortical|telenceph|forebrain", meta$region_clean, ignore.case=TRUE)

bins <- list(
  early_7_8=c(7.0,8.5),
  nine=c(8.5,10.0),
  eleven_twelve=c(10.5,12.5),
  fifteen_seventeen=c(14.5,17.5),
  candidate_11_17=c(11.0,17.5)
)

strict_stats <- function(xy,xx) {
  xy <- xy[is.finite(xy)]; xx <- xx[is.finite(xx)]
  if (!length(xy) || !length(xx)) return(NULL)
  if (min(xy) > max(xx)) return(list(direction="XY>XX", margin=min(xy)-max(xx),
    xy_min=min(xy),xy_max=max(xy),xx_min=min(xx),xx_max=max(xx)))
  if (max(xy) < min(xx)) return(list(direction="XY<XX", margin=min(xx)-max(xy),
    xy_min=min(xy),xy_max=max(xy),xx_min=min(xx),xx_max=max(xx)))
  NULL
}

loo <- function(xy,xx) {
  checks <- c()
  if (length(xy)>2) for (i in seq_along(xy)) checks <- c(checks, !is.null(strict_stats(xy[-i],xx)))
  if (length(xx)>2) for (i in seq_along(xx)) checks <- c(checks, !is.null(strict_stats(xy,xx[-i])))
  if (!length(checks)) return(NA_real_)
  mean(checks)
}

exact_sep_null <- function(nxy,nxx) {
  # If all observations are continuously exchangeable under the null, there are exactly
  # two labelings that put one group entirely above the other.
  2 / choose(nxy+nxx, nxy)
}

donor_matrix <- function(region, lo, hi) {
  region_mask <- if (region == "__CORTEX_BROAD__") meta$is_cortex_broad else meta$region_clean == region
  ii <- which(region_mask & meta$age_wpc >= lo & meta$age_wpc <= hi)
  if (!length(ii)) return(NULL)
  mm <- meta[ii,,drop=FALSE]
  ex <- lcpm[,ii,drop=FALSE]
  key <- paste(mm$donor_id,mm$sex_std,sep="|||")
  # mean within donor+sex when multiple technical/anatomical samples remain under same region label
  keys <- unique(key)
  mat <- sapply(keys, function(k) rowMeans(ex[,key==k,drop=FALSE]))
  if (is.vector(mat)) mat <- matrix(mat,ncol=1,dimnames=list(rownames(ex),keys))
  parts <- strsplit(colnames(mat),"\\|\\|\\|")
  info <- data.frame(donor_id=vapply(parts,`[`,character(1),1),
                     sex=vapply(parts,`[`,character(1),2), stringsAsFactors=FALSE)
  list(expr=mat, info=info)
}

MIN_PER_SEX <- 4L
all_rows <- list()
regions <- c("__CORTEX_BROAD__", sort(unique(meta$region_clean)))
cat("Resolved regions:", paste(regions, collapse=" | "), "\n")

for (region in regions) {
  for (bn in names(bins)) {
    lo <- bins[[bn]][1]; hi <- bins[[bn]][2]
    dm <- donor_matrix(region,lo,hi)
    if (is.null(dm)) next
    nxy <- sum(dm$info$sex=="XY"); nxx <- sum(dm$info$sex=="XX")
    if (nxy < MIN_PER_SEX || nxx < MIN_PER_SEX) next
    for (gi in seq_len(nrow(dm$expr))) {
      xy <- dm$expr[gi,dm$info$sex=="XY"]
      xx <- dm$expr[gi,dm$info$sex=="XX"]
      st <- strict_stats(xy,xx)
      if (is.null(st)) next
      all_rows[[length(all_rows)+1]] <- data.frame(
        region=region,bin=bn,lo_wpc=lo,hi_wpc=hi,gene_id=rownames(dm$expr)[gi],
        n_xy=nxy,n_xx=nxx,direction=st$direction,margin=st$margin,
        xy_min=st$xy_min,xy_max=st$xy_max,xx_min=st$xx_min,xx_max=st$xx_max,
        loo_stability=loo(xy,xx),null_complete_sep_prob=exact_sep_null(nxy,nxx),
        stringsAsFactors=FALSE
      )
    }
  }
}

strict <- if (length(all_rows)) do.call(rbind,all_rows) else data.frame()
write.table(strict, "results/stage1_strict_all.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Temporal emergence: overlap in early_7_8 in the same region, strict later.
temporal_rows <- list()
if (nrow(strict)) {
  later <- strict[strict$bin %in% c("eleven_twelve","fifteen_seventeen","candidate_11_17"),,drop=FALSE]
  for (ri in seq_len(nrow(later))) {
    rr <- later[ri,]
    early <- donor_matrix(rr$region, 7.0, 8.5)
    if (is.null(early)) next
    nxy <- sum(early$info$sex=="XY"); nxx <- sum(early$info$sex=="XX")
    if (nxy < MIN_PER_SEX || nxx < MIN_PER_SEX) next
    g <- rr$gene_id
    if (!g %in% rownames(early$expr)) next
    xy <- early$expr[g,early$info$sex=="XY"]
    xx <- early$expr[g,early$info$sex=="XX"]
    overlap <- !(min(xy) > max(xx) || max(xy) < min(xx))
    if (overlap) {
      rr$early_n_xy <- nxy
      rr$early_n_xx <- nxx
      rr$early_xy_min <- min(xy); rr$early_xy_max <- max(xy)
      rr$early_xx_min <- min(xx); rr$early_xx_max <- max(xx)
      temporal_rows[[length(temporal_rows)+1]] <- rr
    }
  }
}
temporal <- if (length(temporal_rows)) do.call(rbind,temporal_rows) else data.frame()
write.table(temporal, "results/stage1_temporal_emergence.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Add a compact machine-readable summary.
summary <- data.frame(
  metric=c("samples_resolved","donors_resolved","autosomal_genes","strict_rows","temporal_rows"),
  value=c(nrow(meta), length(unique(meta$donor_id)), nrow(lcpm), nrow(strict), nrow(temporal))
)
write.table(summary, "results/stage1_summary.tsv", sep="\t", quote=FALSE, row.names=FALSE)

cat("\nStage 1 complete\n")
print(summary)
if (nrow(temporal)) {
  cat("\nTop temporal candidates:\n")
  print(utils::head(temporal[order(-temporal$loo_stability,-temporal$margin),],20))
}
