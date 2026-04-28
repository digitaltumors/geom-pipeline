suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
  library(yaml)
})

set_vector_memory_limit <- function() {
  if (!exists("mem.maxVSize", mode = "function")) {
    return(invisible(NULL))
  }
  target_gb <- Sys.getenv("PIPELINE_R_MAX_VSIZE_GB", "64")
  target_mb <- suppressWarnings(as.numeric(target_gb) * 1024)
  if (tolower(target_gb) %in% c("inf", "infinite")) {
    target_mb <- Inf
  }
  if (is.na(target_mb)) {
    return(invisible(NULL))
  }
  current_mb <- mem.maxVSize()
  if (
    !is.finite(current_mb) || is.infinite(target_mb) || current_mb < target_mb
  ) {
    try(mem.maxVSize(target_mb), silent = TRUE)
  }
  invisible(NULL)
}

set_vector_memory_limit()

stopf <- function(message, ...) {
  stop(sprintf(message, ...), call. = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(
    path,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    na.strings = c("", "NA"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_csv <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    na.strings = c("", "NA"),
    stringsAsFactors = FALSE
  )
}

normalize_blank_to_na <- function(values) {
  values <- as.character(values)
  blank <- !is.na(values) & !nzchar(trimws(values))
  values[blank] <- NA_character_
  values
}

first_non_missing <- function(...) {
  values <- list(...)
  out <- rep(NA_character_, length(values[[1L]]))
  for (value in values) {
    value <- normalize_blank_to_na(value)
    replace <- is.na(out) & !is.na(value)
    out[replace] <- value[replace]
  }
  out
}

as_logical_flag <- function(values) {
  if (is.logical(values)) {
    return(values)
  }
  lowered <- tolower(trimws(as.character(values)))
  out <- lowered %in% c("true", "t", "1", "yes")
  out[is.na(values) | !nzchar(lowered)] <- FALSE
  out
}

as_integer_or_na <- function(values) {
  suppressWarnings(as.integer(as.character(values)))
}

build_public_drug_metadata <- function(compound_metadata) {
  compound_metadata <- compound_metadata[
    compound_metadata$Source_Subset == "drugs",
    ,
    drop = FALSE
  ]
  if (!nrow(compound_metadata)) {
    stopf("GEOM drug-only metadata table is empty")
  }

  geom_key <- as.character(compound_metadata$Source_SMILES)
  if (any(is.na(geom_key) | !nzchar(trimws(geom_key)))) {
    stopf(
      "GEOM Source_SMILES values must be non-missing before MAE construction"
    )
  }
  if (anyDuplicated(geom_key) > 0L) {
    stopf("GEOM Source_SMILES values must be unique before MAE construction")
  }

  pubchem_cid <- as_integer_or_na(compound_metadata$AnnotationDB_CID)
  public <- data.frame(
    GEOM.Source.SMILES = geom_key,
    Pubchem.CID = pubchem_cid,
    InChIKey = as.character(compound_metadata$Metadata_InChIKey),
    GEOM.RDKit.SMILES = as.character(compound_metadata$Metadata_SMILES),
    In.AnnotationDB = as_logical_flag(compound_metadata$In_AnnotationDB) |
      !is.na(pubchem_cid),
    AnnotationDB.Name = as.character(compound_metadata$AnnotationDB_Name),
    AnnotationDB.SMILES = as.character(compound_metadata$AnnotationDB_SMILES),
    GEOM.Source.Subset = as.character(compound_metadata$Source_Subset),
    GEOM.Canonical.SMILES = as.character(compound_metadata$Canonical_SMILES),
    GEOM.Charge = suppressWarnings(as.numeric(compound_metadata$Charge)),
    GEOM.Total.Conformers = as_integer_or_na(
      compound_metadata$Total_Conformers
    ),
    GEOM.Unique.Conformers = as_integer_or_na(
      compound_metadata$Unique_Conformers
    ),
    GEOM.Lowest.Energy = suppressWarnings(
      as.numeric(compound_metadata$Lowest_Energy)
    ),
    GEOM.Ensemble.Energy = suppressWarnings(
      as.numeric(compound_metadata$Ensemble_Energy)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  row.names(public) <- public$GEOM.Source.SMILES
  public
}

build_whim_experiment <- function(path, assay_name, geom_keys) {
  frame <- read_csv(path)
  if (!("Compound_ID" %in% names(frame))) {
    stopf("%s is missing Compound_ID", path)
  }
  frame$GEOM.Source.SMILES <- as.character(frame$Compound_ID)
  if (anyDuplicated(frame$GEOM.Source.SMILES) > 0L) {
    stopf("%s contains duplicate Compound_ID values", path)
  }
  frame <- frame[match(geom_keys, frame$GEOM.Source.SMILES), , drop = FALSE]
  if (any(is.na(frame$GEOM.Source.SMILES))) {
    stopf("%s is missing rows for selected GEOM source SMILES", path)
  }

  feature_columns <- setdiff(
    names(frame),
    c("Compound_ID", "GEOM.Source.SMILES")
  )
  assay_matrix <- t(as.matrix(frame[, feature_columns, drop = FALSE]))
  storage.mode(assay_matrix) <- "numeric"
  colnames(assay_matrix) <- frame$GEOM.Source.SMILES
  rownames(assay_matrix) <- feature_columns

  row_data <- data.frame(
    Feature.ID = feature_columns,
    Feature.Index = seq_along(feature_columns),
    Assay.Name = assay_name,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  row.names(row_data) <- row_data$Feature.ID

  SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(values = assay_matrix),
    rowData = S4Vectors::DataFrame(row_data, row.names = row.names(row_data)),
    colData = S4Vectors::DataFrame(row.names = colnames(assay_matrix))
  )
}

compound_metadata <- read_tsv(snakemake@input[["compound_metadata"]])
drug_metadata <- build_public_drug_metadata(compound_metadata)
geom_keys <- drug_metadata$GEOM.Source.SMILES

experiments <- list(
  WHIM_hp = build_whim_experiment(
    snakemake@input[["whim_hp"]],
    "WHIM_hp",
    geom_keys
  ),
  WHIMS_avg = build_whim_experiment(
    snakemake@input[["whims_avg"]],
    "WHIMS_avg",
    geom_keys
  ),
  WHIMS_wavg = build_whim_experiment(
    snakemake@input[["whims_wavg"]],
    "WHIMS_wavg",
    geom_keys
  )
)

sample_map <- do.call(
  rbind,
  lapply(names(experiments), function(assay_name) {
    data.frame(
      assay = assay_name,
      primary = colnames(experiments[[assay_name]]),
      colname = colnames(experiments[[assay_name]]),
      stringsAsFactors = FALSE
    )
  })
)

col_data <- drug_metadata
row.names(col_data) <- col_data$GEOM.Source.SMILES

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = MultiAssayExperiment::ExperimentList(experiments),
  colData = S4Vectors::DataFrame(col_data, row.names = row.names(col_data)),
  sampleMap = S4Vectors::DataFrame(sample_map)
)

metadata(mae) <- list(
  Pipeline = list(
    ID = snakemake@params[["dataset_id"]],
    Version = snakemake@params[["dataset_version"]],
    Config = yaml::read_yaml(
      snakemake@input[["configfile"]],
      eval.expr = FALSE
    )
  ),
  Drug.Metadata = S4Vectors::DataFrame(
    drug_metadata,
    row.names = drug_metadata$GEOM.Source.SMILES
  )
)

output_path <- snakemake@output[["mae_rds"]]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(mae, output_path)
