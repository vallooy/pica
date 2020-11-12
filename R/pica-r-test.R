check <- function(x) {
  x <- as.numeric(x)
  if (x > 0) {
    result <- "Positive"
  }
  else if (x < 0) {
    result <- "Negative"
  }
  else {
    result <- "Zero"
  }
  list(
    message = paste(result)
  )
}
