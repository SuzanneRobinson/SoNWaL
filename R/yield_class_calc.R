
#' yield_class_calc: calculate yield class using MAI and CAI from SoNWaL output
#' @param out SoNWaL output
#' @return yield class
#' @export
yield_class_calc <- function(out) {
  #calculate mean annual incrememnt
  MAI <- (aggregate(out$MAI ~ out$Year, FUN = mean))
  #calculate cumulative annual increment
  CAI <- (aggregate(out$CAI ~ out$Year, FUN = mean))
  names(MAI) <- c("x", "y")
  names(CAI) <- c("x", "y")
  MAIx <- MAI[-c(1:10), ]
  CAIx <- CAI[-c(1:10), ]

  MAI$x <- c(1:nrow(MAI))
  CAI$x <- c(1:nrow(CAI))

  CAI$y <- (predict(loess(CAI$y ~ CAI$x, span = 1)))
  MAI[1, 2] <- 0
  MAI[MAI$y==Inf,]<-0
  MAI$y <- (predict(loess(MAI$y ~ MAI$x)))

  #identify point where MAI and CAI intersect to get YC
  inter_sec <- tryCatch({
    curve_intersect(MAI[-c(1:10), ], CAI[-c(1:10), ],
                    empirical = TRUE, domain = NULL)[[2]]

  },
  error = function(cond) {
    CAIx$MAI <- MAIx$y
    CAIx[which(CAIx$y < CAIx$MAI), ][1, 3]
  })

  inter_sec <- if (length(inter_sec) > 1)
    inter_sec[2]
  else
    inter_sec
  return(inter_sec)

}