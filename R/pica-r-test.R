check <- function(x) {
  x <- as.numeric(x)
  if (x > 0) {
    result <- "Guilhem, ceci est positif!"
  }
  else if (x < 0) {
    result <- "Guilhem, ceci est négatif!"
  }
  else {
    result <- "Guilhem, ça c'est nul!"
  }
  list(
    message = paste(result)
  )
}
