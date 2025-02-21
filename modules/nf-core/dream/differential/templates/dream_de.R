#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE){

    ## Get file extension
    ext <- tools::file_ext(basename(file))

    ## Define separator
    if (ext == "tsv" || ext == "txt") {
        separator <- '\t'
    } else if (ext == "csv") {
        separator <- ','
    } else {
        stop(paste("Unknown separator for", ext))
    }

    ## Read file
    cat("Reading file", basename(file), "with", ext, "separator\n")

    df <- read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )

    return(df)
}

parse_args <- function(x) {
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  args_vals <- lapply(args_list, function(x) scan(text = x, what = 'character', quiet = TRUE))
  args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z })
  parsed_args <- structure(lapply(args_vals, function(x) x[2]),
                           names = lapply(args_vals, function(x) x[1]))
  parsed_args[ !is.na(parsed_args) ]
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################
opt <- list(
  output_prefix      = "$meta.contrast_id", # Prefix for output files
  count_file         = "$counts", # File containing raw counts
  sample_file        = "$samplesheet", # File containing sample information
  contrast_variable  = "$meta.contrast_variable", # Variable for contrast
  reference_level    = "$meta.contrast_reference", # Reference level for the contrast
  target_level       = "$meta.contrast_target", # Target level for the contrast
  blocking_variables = "$meta.blocking_factors", # Blocking variables for the analysis
  sample_id_col      = "experiment_accession", # Column name for sample IDs
  threads            = "$task.cpus", # Number of threads for multithreading
  subset_to_contrast_samples            = FALSE, # Whether to subset to contrast samples
  exclude_samples_col = NULL, # Column for excluding samples
  exclude_samples_values = NULL, # Values for excluding samples
  number             = Inf, # Maximum number of results
  adjust.method      = "BH", # Adjustment method for topTable
  p.value            = 1, # P-value threshold for topTable
  lfc                = 0, # Log fold-change threshold for topTable
  confint            = FALSE, # Whether to compute confidence intervals in topTable
  ndups              = NULL, # Number of duplicates for lmFit
  spacing            = NULL, # Spacing for lmFit
  block              = NULL, # Block design for lmFit
  correlation        = NULL, # Correlation for lmFit
  method             = "ls", # Method for lmFit
  proportion         = 0.01, # Proportion for eBayes
  stdev_coef_lim     = "0.1,4", # Standard deviation coefficient limits for eBayes
  trend              = FALSE, # Whether to use trend in eBayes
  robust             = FALSE, # Whether to use robust method in eBayes
  winsor_tail_p      = "0.05,0.1" # Winsor tail probabilities for eBayes
)

args_opt <- parse_args("$task.ext.args")

for (ao in names(args_opt)) {
  if (!ao %in% names(opt)) {
    stop(paste("Invalid option:", ao))
  }
  opt[[ao]] <- args_opt[[ao]]
}

# Check if required parameters have been provided
required_opts <- c("contrast_variable", "reference_level", "target_level", "output_prefix", "count_file", "sample_file")
missing <- required_opts[sapply(required_opts, function(o) is.null(opt[[o]]))]

if (length(missing) > 0) {
    stop(paste("Missing required arguments:", paste(missing, collapse = ", ")), call. = FALSE)
}

# Check if required file inputs are valid
for (file_input in c("count_file", "sample_file")) {
    if (!file.exists(opt[[file_input]])) {
        stop(paste0("Value of ", file_input, ": ", opt[[file_input]], " is not a valid file"), call. = FALSE)
    }
}

# Convert specific options to numeric vectors
vector_opt <- c("stdev_coef_lim", "winsor_tail_p")
opt[vector_opt] <- lapply(strsplit(unlist(opt[vector_opt]), ","), as.numeric)

# Convert string "null" to NULL
if (!is.null(opt\$blocking_variables) && tolower(opt\$blocking_variables) == "null") {
    opt\$blocking_variables <- NULL
}

# Save first version of RData, useful for debuging with all original parameters already set
save.image("dream_de.RData")

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(edgeR)
library(variancePartition)

## Load RData (for dev purposes)
#setwd("/workspace/differentialabundance/results/work/fe/264e25e728b3300231e7d393e939f4")
#load("dream_de.RData")

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

intensities.table <- read_delim_flexible(file = opt\$count_file)
sample.sheet <- read_delim_flexible(file = opt\$sample_file)

# Deal with spaces that may be in sample column
opt\$sample_id_col <- make.names(opt\$sample_id_col)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

# Check that all samples specified in the input sheet are present in the
# intensities table. Assuming they are, subset and sort the count table to
# match the sample sheet

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(intensities.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else {
    # Save any non-count data, will gene metadata etc we might need later
    nonintensities.table <-
        intensities.table[, !colnames(intensities.table) %in% rownames(sample.sheet), drop = FALSE]
    intensities.table <- intensities.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

contrast_variable <- make.names(opt\$contrast_variable)
blocking.vars <- c()

if (!contrast_variable %in% colnames(sample.sheet)) {
    stop(
        paste0(
        'Chosen contrast variable "',
        contrast_variable,
        '" not in sample sheet'
        )
    )
} else if (any(!c(opt\$reference_level, opt\$target_level) %in% sample.sheet[[contrast_variable]])) {
    stop(
        paste0(
        'Please choose reference and target levels that are present in the ',
        contrast_variable,
        ' column of the sample sheet'
        )
    )
} else if (!is.null(opt\$blocking_variables)) {
    blocking.vars = make.names(unlist(strsplit(opt\$blocking_variables, split = ';')))
    if (!all(blocking.vars %in% colnames(sample.sheet))) {
        missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
        stop(
            paste(
                'Blocking variables', missing_block,
                'do not correspond to sample sheet columns.'
            )
        )
    }
}

# Handle conflicts between blocking variables and block
if (!is.null(opt\$block) && !is.null(opt\$blocking_variables)) {
    if (opt\$block %in% blocking.vars) {
        warning(paste("Variable", opt\$block, "is specified both as a fixed effect and a random effect. It will be treated as a random effect only."))
        blocking.vars <- setdiff(blocking.vars, opt\$block)
        if (length(blocking.vars) == 0) {
            opt\$blocking_variables <- NULL
        } else {
            opt\$blocking_variables <- paste(blocking.vars, collapse = ';')
        }
    }
}

# Optionally, subset to only the samples involved in the contrast

if (opt\$subset_to_contrast_samples){
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt\$target_level, opt\$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    intensities.table <- intensities.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((! is.null(opt\$exclude_samples_col)) && (! is.null(opt\$exclude_samples_values))){
    exclude_values = unlist(strsplit(opt\$exclude_samples_values, split = ';'))

    if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
    sample_selector <- ! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values

    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    intensities.table <- intensities.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

################################################
################################################
## Build the Model Formula                    ##
################################################
################################################

# Build the model formula with blocking variables first
model_vars <- c()

if (!is.null(opt\$blocking_variables)) {
    # Include blocking variables (including pairing variables if any)
    model_vars <- c(model_vars, blocking.vars)
}

# Add the contrast variable at the end
model_vars <- c(model_vars, contrast_variable)

# Construct the model formula
model <- paste('~ 0 +', paste(model_vars, collapse = '+'))

# Make sure all the appropriate variables are factors
vars_to_factor <- model_vars  # All variables in the model need to be factors
for (v in vars_to_factor) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
}

################################################
################################################
## Run Dream processes                        ##
################################################
################################################

# Generate the design matrix
design <- model.matrix(
    as.formula(model),
    data=sample.sheet
)


# Specify parallel processing
param <- SnowParam(as.numeric(opt\$threads), "SOCK", progressbar = TRUE)

## Set formula
#form <- ~ Disease + (1 | Individual)
form <- as.formula(model)

# Create a DGEList object for RNA-seq data
dge <- DGEList(counts = intensities.table)

## Calculate normalization factors
dge <- calcNormFactors(dge)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, sample.sheet, BPPARAM = param)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation

## TODO: this is the placeholder for `L = makeContrastsDream(...)` function, whose object should be included with `L` argument in dream() function
L <- variancePartition::makeContrastsDream(form, sample.sheet,
    contrasts = c(
        setNames(paste(
            paste0(opt\$contrast_variable, opt\$target_level),
            paste0(opt\$contrast_variable, opt\$reference_level),
            sep = " - "), opt\$output_prefix)
        )
    )

# Visualize contrast matrix
contrasts_plot <- plotContrasts(L)

png(
    file = paste(opt\$output_prefix, 'dream.contrasts_plot.png', sep = '.'),
    width = 600,
    height = 300
)
plot(contrasts_plot)
dev.off()

# fit dream model with contrasts
fitmm <- dream(vobjDream, form, sample.sheet, L)
fitmm <- variancePartition::eBayes(fitmm)

# get names of available coefficients and contrasts for testing
colnames(fitmm)

# Get results of hypothesis test on coefficients of interest
for (COEFFICIENT in opt\$output_prefix) {

    ## Initialize topTable() arguments
    toptable_args <- list(
        fit = fitmm,
        coef = COEFFICIENT,
        sort.by = 'none',
        number = nrow(intensities.table)
    )

    ## Complete list with extra arguments if they were provided
    if (! is.null(opt\$adjust.method)){
        toptable_args[['adjust.method']] <- opt\$adjust.method
    }
    if (! is.null(opt\$p.value)){
        toptable_args[['p.value']] <- as.numeric(opt\$p.value)
    }
    if (! is.null(opt\$lfc)){
        toptable_args[['lfc']] <- as.numeric(opt\$lfc)
    }
    if (! is.null(opt\$confint)){
        toptable_args[['confint']] <- as.logical(opt\$confint)
    }

    ## generate topTable
    comp.results <- do.call(topTable, toptable_args)[rownames(intensities.table),]

    ## Export topTable
    write.table(
        comp.results,
        file = paste(opt\$output_prefix, 'dream.results.tsv', sep = '.'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

## END OF NEW

################################################
################################################
## Generate outputs                           ##
################################################
################################################

# Dispersion plot

png(
    file = paste(opt\$output_prefix, 'dream.mean_difference.png', sep = '.'),
    width = 600,
    height = 600
)
plotMD(fitmm)
dev.off()

# R object for other processes to use
saveRDS(fitmm, file = paste(opt\$output_prefix, 'MArrayMM.dream.rds', sep = '.'))

# Save model to file
write(model, file=paste(opt\$output_prefix, 'dream.model.txt', sep = '.'))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt\$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
################################################
################################################

save.image("dream_de.RData")

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
edger.version <- as.character(packageVersion('edgeR'))
variancePartition.version <- as.character(packageVersion('variancePartition'))


writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-edger:', edger.version),
        paste('    variancePartition:', variancePartition.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
