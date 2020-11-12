check <- function(x) {
  x <- as.numeric(x)
  if (x > 0) {
    result <- "Guilhem, selon R, ceci est positif!"
  }
  else if (x < 0) {
    result <- "Guilhem, selon R, ceci est negatif!"
  }
  else {
    result <- "Guilhem, selon R, c'est nul!"
  }
  list(message = paste(result))
}
