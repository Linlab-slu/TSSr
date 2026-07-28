canonicalize_tssr_value <- function(value) {
    if (data.table::is.data.table(value)) {
        return(list(
            data = as.data.frame(data.table::copy(value)),
            key = data.table::key(value)
        ))
    }
    if (is.list(value)) {
        return(lapply(value, canonicalize_tssr_value))
    }
    value
}

tssr_content <- function(object) {
    slots <- methods::slotNames(object)
    stats::setNames(
        lapply(slots, function(slot_name) {
            canonicalize_tssr_value(methods::slot(object, slot_name))
        }),
        slots
    )
}

expect_tssr_content_equal <- function(object, expected) {
    testthat::expect_equal(tssr_content(object), expected)
}
