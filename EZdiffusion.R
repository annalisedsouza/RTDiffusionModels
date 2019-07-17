# Code obtained from Wagenmakers, van der Maas, and Grasman (2007)

get.vaTer <- function(Pc, VRT, MRT, s=0.1)
  {
  s2 <- s^2
  # The default value for the scaling parameter s equals 0.1
  if (Pc == 0)
    cat("Oops, Pc == 0!\n")
  if (Pc == 1)
    cat("Oops, Pc == 1!\n")
  # If Pc equals 0, .5, or 1, the method will not work, and
  # an edge correction is required.
  L <- qlogis(Pc)
  # The function “qlogis” calculates the logit.
  x <- L*(L*Pc^2 - L*Pc + Pc - .5)/VRT
  v <- sign(Pc-.5)*s*x^(1/4)
  # This gives drift rate.
  a <- s2*qlogis(Pc)/v
  # This gives boundary separation.
  y <- -v*a/s2
  MDT <- (a/(2*v)) * (1-exp(y))/(1+exp(y))
  Ter <- MRT - MDT
  # This gives nondecision time.
  return(list(v, a, Ter))
}
