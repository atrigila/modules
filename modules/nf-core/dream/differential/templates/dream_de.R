#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################
cat("Initializing functions\n")

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



isMixedModelFormula <- function(formula) {

    !is.null(lme4::findbars(as.formula(formula)))

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
  winsor_tail_p      = "0.05,0.1", # Winsor tail probabilities for eBayes
  ddf                = "adaptive", # "Specifiy 'Satterthwaite', 'Kenward-Roger', or 'adaptive' method for dream()
  reml               = FALSE
)

args_opt <- parse_args("$task.ext.args")

for (ao in names(args_opt)) {
  if (!ao %in% names(opt)) {
    stop(paste("Invalid option:", ao))
  }
  opt[[ao]] <- args_opt[[ao]]
}

# Check if required parameters have been provided
cat("Validating required arguments\n")

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

## Check default values for ddf options
ddf_valid <- c("Satterthwaite", "Kenward-Roger", "adaptive")
if ( !opt\$ddf %in% ddf_valid ) {
    stop(paste0("'--ddf '", opt\$ddf, "' is not a valid option from '", paste(ddf_valid, collapse = "', '"), "'"), call. = FALSE)
}

## Check adjust method options
adjust_valid <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
if (!opt\$adjust.method %in% adjust_valid && !is.null(opt\$adjust.method)) {
        stop(paste0("'--adjust.method '", opt\$adjust.method, "' is not a valid option from '", paste(adjust_valid, collapse = "', '"), "', or NULL"), call. = FALSE)
}

# Convert specific options to numeric vectors
vector_opt <- c("stdev_coef_lim", "winsor_tail_p")
opt[vector_opt] <- lapply(strsplit(unlist(opt[vector_opt]), ","), as.numeric)

# Convert string "null" to NULL
if (!is.null(opt\$blocking_variables) && tolower(opt\$blocking_variables) == "null") {
    opt\$blocking_variables <- NULL
}

# Save first version of RData, useful for debuging with all original parameters already set
cat("Exporting preliminary RData\n")

work_dir <- getwd()                         ## for dev purposes
save.image("dream_de.RData")
#setwd(work_dir); load("dream_de.RData")    ## for dev purposes

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################
cat("Importing libraries\n")

library(edgeR)
library(variancePartition)

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
cat("Validating contrasts\n")

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

## TODO: START ADAPTING TO DREAM IN HERE!
## NEW STARTS HERE!
cat("Creating formula\n")

# Build the model formula with blocking variables first
model_vars <- c()

if (!is.null(opt\$blocking_variables)) {
    cat("opt\$blocking_variables:", paste(opt\$blocking_variables, collapse = ' '), "\n")
    print(opt\$blocking_variables)

    # Include blocking variables (including pairing variables if any)
    for (VARIABLE in opt\$blocking_variables) {
        cat("Adding", VARIABLE, "factor to formula\n")

        # Convert to factor
        sample.sheet[[VARIABLE]] <- as.factor(sample.sheet[[VARIABLE]])

        ## Create a string to reconstruct Wilkinson formula " (1 | variable )"
        model_vars <- c(model_vars, paste0("(1 | ", VARIABLE, ")"))

        ## Convert the variable into factor
        if (!is.numeric(sample.sheet[[ VARIABLE ]])) {
            cat("Converting into factor\n")
            sample.sheet[[ VARIABLE ]] <- as.factor(sample.sheet[[ VARIABLE ]])
        }
    }
}

# Construct the model formula
## Expected structure:
## "~ 0 + fixed_effect + (1 | random_variable_1) + (1 | random_variable_N)"
model <- paste(
    '~ 0 +',
    contrast_variable,
    if (length(model_vars) > 0) { paste0(" +", paste(model_vars, collapse = ' + ')) } else { NULL },
    sep = "")  ## TODO: This is limited to additive models! not possible for interaction relations

cat("Model vars:", model_vars, "\n")
cat("Model:", model, "\n")

# Construct the formula
form <- as.formula(model)
cat("Formula:", deparse(form), "\n")

## TODO: Check if the model is mixed or not, could be useful to report it later
cat("Checking for mixed formula\n")
mixed_form <- isMixedModelFormula(form)

################################################
################################################
## Run Dream processes                        ##
################################################
################################################

# Generate the design matrix
sample.sheet[[contrast_variable]] <- as.factor(sample.sheet[[contrast_variable]])
cat("Creating design matrix\n")

design <- model.matrix(
    form,
    sample.sheet
)

# Specify parallel processing
param <- SnowParam(as.numeric(opt\$threads), "SOCK", progressbar = TRUE)

# Create a DGEList object for RNA-seq data
cat("Creating DGEList\n")
dge <- DGEList(counts = intensities.table)

## Calculate normalization factors
cat("Calculating normalization factors\n")
dge <- calcNormFactors(dge)

# estimate weights using linear mixed model of dream
cat("Normalizing data\n")
vobjDream <- voomWithDreamWeights(dge, form, sample.sheet, BPPARAM = param)

# Create and export variance plot
cat("Analyzing variance\n")
if (!is.null(opt\$blocking_variables)) {
  # When blocking factors exist, remove treatment so that all categorical variables are random
  vp_formula <- update(form, ~ . - contrast_variable + 1)
} else {
  # When no blocking factors exist, use the original formula, but add intercept because it is required
  vp_formula <- update(form, ~ . + 1)
}

vp <- fitExtractVarPartModel(
  exprObj = vobjDream,
  formula = vp_formula,
  data = sample.sheet,
  REML = opt\$reml
)

cat("Creating variance plot\n")
var_plot <- plotVarPart(sortCols(vp))

cat("Exporting variance plot\n")
png(
    file = paste(opt\$output_prefix, 'dream.var_plot.png', sep = '.'),
    width = 600,
    height = 300
)
plot(var_plot)
dev.off()

## Set contrast (this can be scaled for more than one comparison)
cat("Building contrasts\n")
L <- variancePartition::makeContrastsDream(
    form,
    sample.sheet,
    contrasts = c(                                                      ## This is a named vector for all the contrast we'd like to run.
        setNames(                                                       ## For now, it's just one contrast but we should be able to scale it with the yml somehow
            paste(                                                      ## For a variable named `treatment`, with levels A and B
                paste0(opt\$contrast_variable, opt\$target_level),            ## Create target value: treatmentB
                paste0(opt\$contrast_variable, opt\$reference_level),         ## Create reference value: treatmentA
                sep = " - "),                                               ## Set a difference between them
            opt\$output_prefix)                                          ## Assign the name to the comparison (`contrast id`` from the yml)
        )
    )

# Visualize contrast matrix
cat("Creating and exporting contrast plot\n")
contrasts_plot <- plotContrasts(L)

# Export contrast plot
png(
    file = paste(opt\$output_prefix, 'dream.contrasts_plot.png', sep = '.'),
    width = 600,
    height = 300
)
plot(contrasts_plot)
dev.off()


# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation (this is the behavior indicated with `ddf = adaptive`)
# We can force it to use the others methods with is ddf argument
cat("Fitting model with dream()\n")
fitmm <-
    dream(
        exprObj = vobjDream,
        formula = form,
        data    = sample.sheet,
        L       = L,
        ddf     = opt\$ddf,
        reml    = opt\$reml
    )

# Adjust results with Empirical Bayes
## Create list of arguments for eBayes
cat("Adjusting results with eBayes\n")
ebayes_args <- list(
    fit = fitmm
)

if (! is.null(opt\$proportion)){
    ebayes_args[['proportion']] <- as.numeric(opt\$proportion)
}
if (! is.null(opt\$stdev.coef.lim)){
    ebayes_args[['stdev.coef.lim']] <- as.numeric(opt\$stdev_coef_lim)
}
if (! is.null(opt\$trend)){
    ebayes_args[['trend']] <- as.logical(opt\$trend)
}
if (! is.null(opt\$robust)){
    ebayes_args[['robust']] <- as.logical(opt\$robust)
}
if (! is.null(opt\$winsor.tail.p)){
    ebayes_args[['winsor.tail.p']] <- as.numeric(opt\$winsor_tail_p)
}

## Run variancePartition::eBayes
fitmm <- do.call(variancePartition::eBayes, ebayes_args)

# get names of available coefficients and contrasts for testing
colnames(fitmm)

# Get results of hypothesis test on coefficients of interest (only one coeff for now)
cat("Exporting results with topTable()\n")
for (COEFFICIENT in opt\$output_prefix) {

    ## Initialize topTable() arguments
    toptable_args <- list(
        fit = fitmm,
        coef = COEFFICIENT,
        sort.by = 'none',
        number = nrow(intensities.table)
    )

    ## Complete list with extra arguments, if they were provided
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
    comp.results <- do.call(variancePartition::topTable, toptable_args)[rownames(intensities.table),]

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
## Generate other outputs                     ##
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
