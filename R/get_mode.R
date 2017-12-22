get_mode <- function(v) {
    if (all(length(unique(v)) == 1)) {
        return(unique(v))
    } else {
        d <- density(v)
        d$x[d$y == max(d$y)] %>% .[ceiling(length(.) / 2)]
    }
}
