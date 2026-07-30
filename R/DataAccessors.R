.asAccessorDataFrame <- function(value) {
    value <- data.table::copy(value)
    data.table::setDF(value)
    value
}

.accessNamedTables <- function(tables, name = NULL, name_label = "sample") {
    if (is.null(name)) {
        return(lapply(tables, .asAccessorDataFrame))
    }
    if (!is.character(name) || length(name) != 1L || is.na(name)) {
        stop(name_label, " must be NULL or one non-missing character value.")
    }
    if (is.null(names(tables)) || !name %in% names(tables)) {
        available <- if (length(tables) == 0L) {
            "none"
        } else {
            paste(names(tables), collapse = ", ")
        }
        stop(
            "Unknown ", name_label, " '", name,
            "'. Available values: ", available, "."
        )
    }
    .asAccessorDataFrame(tables[[name]])
}

.accessResultTables <- function(object, slot_name, name = NULL,
                                name_label = "sample") {
    .accessNamedTables(
        methods::slot(object, slot_name),
        name,
        name_label
    )
}

#' Access data stored in a TSSr object
#'
#' These read-only accessors provide stable alternatives to direct slot
#' access. Tables are returned as independent base \code{data.frame} copies,
#' so modifying a returned value does not modify \code{object} or expose its
#' internal \code{data.table} representation.
#'
#' @param object A TSSr object.
#' @param data Which TSS matrix to return: \code{"raw"} or
#'   \code{"processed"}.
#' @param sample A sample name, or \code{NULL} to return every sample.
#' @param comparison A comparison name, or \code{NULL} to return every
#'   comparison.
#' @param result For \code{DEtables()}, return either \code{"all"} tested
#'   genes or only \code{"significant"} genes.
#'
#' @return \code{TSSmatrix()} and \code{refTable()} return a base
#'   \code{data.frame}. \code{librarySizes()} returns a named numeric vector.
#'   The analysis-result accessors return a base \code{data.frame} when a
#'   sample or comparison is selected, and otherwise return a named list of
#'   base \code{data.frame} objects.
#'
#' @name TSSr-accessors
#'
#' @examples
#' data(exampleTSSr)
#'
#' head(TSSmatrix(exampleTSSr, data = "raw"))
#' librarySizes(exampleTSSr)
#' head(refTable(exampleTSSr))
#'
#' head(tagClusters(exampleTSSr, sample = "control"))
#' head(consensusClusters(exampleTSSr, sample = "control"))
#' head(clusterShape(exampleTSSr, sample = "control"))
#' head(assignedClusters(exampleTSSr, sample = "control"))
#' head(unassignedClusters(exampleTSSr, sample = "control"))
#' head(filteredClusters(exampleTSSr, sample = "control"))
#' head(enhancers(exampleTSSr, sample = "control"))
#' head(DEtables(
#'     exampleTSSr,
#'     comparison = "control_VS_treat",
#'     result = "significant"
#' ))
#' head(TAGtables(exampleTSSr, sample = "control"))
#' head(PromoterShift(exampleTSSr, comparison = "control_VS_treat"))
NULL

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "TSSmatrix",
    function(object, data = c("raw", "processed"))
        standardGeneric("TSSmatrix"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "TSSmatrix",
    signature(object = "TSSr"),
    function(object, data = c("raw", "processed")) {
        data <- match.arg(data)
        value <- if (identical(data, "raw")) {
            object@TSSrawMatrix
        } else {
            object@TSSprocessedMatrix
        }
        .asAccessorDataFrame(value)
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "librarySizes",
    function(object) standardGeneric("librarySizes"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "librarySizes",
    signature(object = "TSSr"),
    function(object) {
        structure(
            as.numeric(object@librarySizes),
            names = names(object@librarySizes)
        )
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "refTable",
    function(object) standardGeneric("refTable"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "refTable",
    signature(object = "TSSr"),
    function(object) .asAccessorDataFrame(object@refTable)
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "tagClusters",
    function(object, sample = NULL) standardGeneric("tagClusters"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "tagClusters",
    signature(object = "TSSr"),
    function(object, sample = NULL) {
        .accessResultTables(object, "tagClusters", sample, "sample")
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "consensusClusters",
    function(object, sample = NULL) standardGeneric("consensusClusters"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "consensusClusters",
    signature(object = "TSSr"),
    function(object, sample = NULL) {
        .accessResultTables(object, "consensusClusters", sample, "sample")
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "clusterShape",
    function(object, sample = NULL) standardGeneric("clusterShape"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "clusterShape",
    signature(object = "TSSr"),
    function(object, sample = NULL) {
        .accessResultTables(object, "clusterShape", sample, "sample")
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "assignedClusters",
    function(object, sample = NULL) standardGeneric("assignedClusters"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "assignedClusters",
    signature(object = "TSSr"),
    function(object, sample = NULL) {
        .accessResultTables(object, "assignedClusters", sample, "sample")
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "unassignedClusters",
    function(object, sample = NULL) standardGeneric("unassignedClusters"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "unassignedClusters",
    signature(object = "TSSr"),
    function(object, sample = NULL) {
        .accessResultTables(object, "unassignedClusters", sample, "sample")
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "filteredClusters",
    function(object, sample = NULL) standardGeneric("filteredClusters"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "filteredClusters",
    signature(object = "TSSr"),
    function(object, sample = NULL) {
        .accessResultTables(object, "filteredClusters", sample, "sample")
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "enhancers",
    function(object, sample = NULL) standardGeneric("enhancers"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "enhancers",
    signature(object = "TSSr"),
    function(object, sample = NULL) {
        .accessResultTables(object, "enhancers", sample, "sample")
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "DEtables",
    function(object, comparison = NULL,
             result = c("all", "significant")) standardGeneric("DEtables"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "DEtables",
    signature(object = "TSSr"),
    function(object, comparison = NULL,
             result = c("all", "significant")) {
        result <- match.arg(result)
        result_name <- if (identical(result, "all")) "DEtable" else "DEsig"
        tables <- lapply(object@DEtables, function(value) {
            value[[result_name]]
        })
        .accessNamedTables(
            tables,
            comparison,
            "comparison"
        )
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "TAGtables",
    function(object, sample = NULL) standardGeneric("TAGtables"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "TAGtables",
    signature(object = "TSSr"),
    function(object, sample = NULL) {
        .accessResultTables(object, "TAGtables", sample, "sample")
    }
)

#' @rdname TSSr-accessors
#' @export
setGeneric(
    "PromoterShift",
    function(object, comparison = NULL) standardGeneric("PromoterShift"),
    signature = "object"
)

#' @rdname TSSr-accessors
#' @export
setMethod(
    "PromoterShift",
    signature(object = "TSSr"),
    function(object, comparison = NULL) {
        .accessResultTables(
            object,
            "PromoterShift",
            comparison,
            "comparison"
        )
    }
)
