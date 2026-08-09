.mapPointsToIntervals <- function(
  points, intervals, pointColumn = "pos",
  intervalStart = "start", intervalEnd = "end"
) {
    chr <- interval_id <- point_id <- NULL

    if (nrow(points) == 0L || nrow(intervals) == 0L) {
        return(data.table(
            interval_id = integer(),
            point_id = integer()
        ))
    }

    interval.work <- data.table(
        chr = as.character(intervals[["chr"]]),
        strand = as.character(intervals[["strand"]]),
        start = intervals[[intervalStart]],
        end = intervals[[intervalEnd]],
        interval_id = seq_len(nrow(intervals))
    )
    point.position <- points[[pointColumn]]
    point.work <- data.table(
        chr = as.character(points[["chr"]]),
        strand = as.character(points[["strand"]]),
        start = point.position,
        end = point.position,
        point_id = seq_len(nrow(points))
    )

    setkey(interval.work, chr, strand, start, end)
    hits <- foverlaps(
        point.work,
        interval.work,
        by.x = c("chr", "strand", "start", "end"),
        by.y = c("chr", "strand", "start", "end"),
        type = "within",
        nomatch = NULL
    )
    hits <- hits[, list(interval_id, point_id)]
    setorder(hits, interval_id, point_id)
    hits
}
