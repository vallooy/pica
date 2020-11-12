check <- function(x) {
  x <- as.numeric(x)
  if (x > 0) {
    result <- "Guilhem, selon R, ceci est positif!"
  }
  else if (x < 0) {
    result <- "Guilhem, selon R, ceci est négatif!"
  }
  else {
    result <- "Guilhem, selon R, ça c'est nul!"
  }
  list(message = paste(result))
}
