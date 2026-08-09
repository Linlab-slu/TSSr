#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 5L) {
    stop(
        paste(
            "Usage: Rscript benchmark-pipeline-stages.R CONFIG OUTDIR",
            "[STAGES] [REPETITIONS] [BASELINE_DIR]"
        ),
        call. = FALSE
    )
}

config.file <- normalizePath(args[[1L]], mustWork = TRUE)
output.dir <- args[[2L]]
stage.argument <- if (length(args) >= 3L) args[[3L]] else "all"
repetitions <- if (length(args) >= 4L) as.integer(args[[4L]]) else 1L
baseline.dir <- if (length(args) >= 5L) {
    normalizePath(args[[5L]], mustWork = TRUE)
} else {
    NULL
}

if (is.na(repetitions) || repetitions < 1L) {
    stop("REPETITIONS must be a positive integer.", call. = FALSE)
}

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
output.dir <- normalizePath(output.dir, mustWork = TRUE)

config.environment <- new.env(parent = globalenv())
sys.source(config.file, envir = config.environment)
config.values <- as.list(config.environment, all.names = TRUE)
invisible(list2env(config.values, envir = environment()))

suppressPackageStartupMessages(library(TSSr))

required.values <- c(
    "input_type", "input_files", "gff_file", "genome_name",
    "sample_labels", "sample_labels_merged", "merge_index",
    "organism_name", "sequencing_quality_threshold",
    "mapping_quality_threshold", "softclipping_allowed",
    "filter_method", "filter_tpm_low", "cluster_method",
    "peak_distance", "extension_distance", "local_threshold",
    "cluster_threshold", "consensus_distance", "shape_method",
    "annotation_type", "filter_cluster", "filter_cluster_threshold",
    "annotation_upstream", "annotation_upstream_overlap",
    "annotation_downstream", "enhancer_flanking",
    "enhancer_distance_to_gene", "comparison_pairs", "de_pvalue",
    "shift_collection_pvalue", "use_multicore", "number_of_cores"
)
missing.values <- setdiff(required.values, names(config.values))
if (length(missing.values) > 0L) {
    stop(
        "Configuration is missing: ",
        paste(missing.values, collapse = ", "),
        call. = FALSE
    )
}

stage.names <- c(
    "getTSS", "merge", "normalize", "filter", "cluster", "consensus",
    "shape", "annotate", "enhancer", "de", "shift"
)
selected.stages <- if (identical(stage.argument, "all")) {
    stage.names
} else {
    trimws(strsplit(stage.argument, ",", fixed = TRUE)[[1L]])
}
unknown.stages <- setdiff(selected.stages, stage.names)
if (length(unknown.stages) > 0L) {
    stop(
        "Unknown stage(s): ", paste(unknown.stages, collapse = ", "),
        call. = FALSE
    )
}
if (is.null(baseline.dir) && !identical(selected.stages, stage.names)) {
    stop(
        "BASELINE_DIR is required when running selected stages.",
        call. = FALSE
    )
}

stage.files <- setNames(
    sprintf("%02d_%s.rds", seq_along(stage.names), stage.names),
    stage.names
)
predecessor.files <- c(
    getTSS = "00_constructed.rds",
    setNames(stage.files[-length(stage.files)], stage.names[-1L])
)

restore.tables <- function(value) {
    if (is.data.frame(value)) {
        return(data.table::as.data.table(value))
    }
    if (is.list(value)) {
        return(lapply(value, restore.tables))
    }
    value
}

read.benchmark.object <- function(path) {
    value <- readRDS(path)
    if (methods::is(value, "TSSr")) {
        return(value)
    }
    if (!is.list(value) || !is.list(value$slots)) {
        stop(
            "Snapshot is neither a TSSr object nor a canonical TSSr record: ",
            path,
            call. = FALSE
        )
    }

    object <- methods::new("TSSr")
    available.slots <- intersect(
        methods::slotNames(object), names(value$slots)
    )
    for (slot.name in available.slots) {
        methods::slot(object, slot.name, check = FALSE) <-
            restore.tables(value$slots[[slot.name]])
    }
    methods::validObject(object)
    object
}

new.object <- function() {
    methods::new(
        "TSSr",
        genomeName = genome_name,
        inputFiles = input_files,
        inputFilesType = input_type,
        sampleLabels = sample_labels,
        sampleLabelsMerged = sample_labels_merged,
        mergeIndex = merge_index,
        refSource = gff_file,
        organismName = organism_name
    )
}

run.operation <- function(stage, object) {
    switch(
        stage,
        getTSS = getTSS(
            object,
            sequencingQualityThreshold = sequencing_quality_threshold,
            mappingQualityThreshold = mapping_quality_threshold,
            softclippingAllowed = softclipping_allowed
        ),
        merge = mergeSamples(object, mergeIndex = merge_index),
        normalize = normalizeTSS(object),
        filter = filterTSS(
            object, method = filter_method, tpmLow = filter_tpm_low
        ),
        cluster = clusterTSS(
            object,
            method = cluster_method,
            peakDistance = peak_distance,
            extensionDistance = extension_distance,
            localThreshold = local_threshold,
            clusterThreshold = cluster_threshold,
            useMultiCore = use_multicore,
            numCores = number_of_cores
        ),
        consensus = consensusCluster(
            object,
            dis = consensus_distance,
            useMultiCore = use_multicore,
            numCores = number_of_cores
        ),
        shape = shapeCluster(
            object,
            clusters = "consensusClusters",
            method = shape_method,
            useMultiCore = use_multicore,
            numCores = number_of_cores
        ),
        annotate = annotateCluster(
            object,
            clusters = "consensusClusters",
            filterCluster = filter_cluster,
            filterClusterThreshold = filter_cluster_threshold,
            annotationType = annotation_type,
            upstream = annotation_upstream,
            upstreamOverlap = annotation_upstream_overlap,
            downstream = annotation_downstream
        ),
        enhancer = callEnhancer(
            object,
            flanking = enhancer_flanking,
            dis2gene = enhancer_distance_to_gene
        ),
        de = deGene(
            object,
            comparePairs = comparison_pairs,
            pval = de_pvalue,
            useMultiCore = use_multicore,
            numCores = number_of_cores
        ),
        shift = shiftPromoter(
            object,
            comparePairs = comparison_pairs,
            pval = shift_collection_pvalue
        )
    )
}

summarize.object <- function(object, stage, iteration, elapsed) {
    data.frame(
        stage = stage,
        iteration = iteration,
        elapsed_seconds = elapsed,
        raw_rows = nrow(object@TSSrawMatrix),
        processed_rows = nrow(object@TSSprocessedMatrix),
        tag_cluster_rows = sum(vapply(
            object@tagClusters, nrow, integer(1)
        )),
        consensus_cluster_rows = sum(vapply(
            object@consensusClusters, nrow, integer(1)
        )),
        stringsAsFactors = FALSE
    )
}

warnings.by.run <- list()
timings <- list()

write.progress <- function() {
    if (length(timings) == 0L) {
        return(invisible(NULL))
    }
    utils::write.table(
        do.call(rbind, timings),
        file.path(output.dir, "timings.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    saveRDS(
        warnings.by.run,
        file.path(output.dir, "warnings.rds"),
        version = 3
    )
    invisible(NULL)
}

if (is.null(baseline.dir)) {
    object <- new.object()
    saveRDS(object, file.path(output.dir, "00_constructed.rds"), version = 3)
    for (stage in stage.names) {
        stage.warnings <- character()
        gc()
        started <- proc.time()[["elapsed"]]
        object <- withCallingHandlers(
            run.operation(stage, object),
            warning = function(condition) {
                stage.warnings <<- c(
                    stage.warnings, conditionMessage(condition)
                )
            }
        )
        elapsed <- proc.time()[["elapsed"]] - started
        saveRDS(object, file.path(output.dir, stage.files[[stage]]), version = 3)
        timings[[length(timings) + 1L]] <- summarize.object(
            object, stage, 1L, elapsed
        )
        warnings.by.run[[paste(stage, 1L, sep = ":")]] <- stage.warnings
        write.progress()
    }
} else {
    for (stage in selected.stages) {
        input.file <- file.path(baseline.dir, predecessor.files[[stage]])
        if (!file.exists(input.file)) {
            stop("Missing predecessor snapshot: ", input.file, call. = FALSE)
        }
        for (iteration in seq_len(repetitions)) {
            object <- read.benchmark.object(input.file)
            stage.warnings <- character()
            gc()
            started <- proc.time()[["elapsed"]]
            result <- withCallingHandlers(
                run.operation(stage, object),
                warning = function(condition) {
                    stage.warnings <<- c(
                        stage.warnings, conditionMessage(condition)
                    )
                }
            )
            elapsed <- proc.time()[["elapsed"]] - started
            timings[[length(timings) + 1L]] <- summarize.object(
                result, stage, iteration, elapsed
            )
            warnings.by.run[[paste(stage, iteration, sep = ":")]] <-
                stage.warnings
            if (iteration == 1L) {
                saveRDS(
                    result,
                    file.path(output.dir, stage.files[[stage]]),
                    version = 3
                )
            }
            write.progress()
        }
    }
}

timings <- do.call(rbind, timings)

metadata <- list(
    git_commit = Sys.getenv("TSSR_BENCHMARK_COMMIT", unset = NA_character_),
    package_version = as.character(packageVersion("TSSr")),
    package_path = find.package("TSSr"),
    config_file = config.file,
    config_md5 = unname(tools::md5sum(config.file)),
    input_files = ifelse(
        file.exists(input_files),
        normalizePath(input_files, mustWork = FALSE),
        input_files
    ),
    input_md5 = ifelse(
        file.exists(input_files),
        unname(tools::md5sum(input_files)),
        NA_character_
    ),
    selected_stages = selected.stages,
    repetitions = repetitions,
    baseline_dir = baseline.dir,
    system_info = Sys.info(),
    session_info = capture.output(sessionInfo()),
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(metadata, file.path(output.dir, "metadata.rds"), version = 3)
