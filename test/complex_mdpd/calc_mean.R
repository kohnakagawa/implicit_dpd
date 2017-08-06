LoadPressure <- function(root.dir) {
  fname <- file.path(root.dir, "pressure.txt")
  pxyz  <- data.matrix(read.table(fname))
  return (pxyz)
}

Slice <- function(vec, begin) {
  return (vec[seq(begin, length(vec))])
}

CalcPressure <- function(pxyz, begin) {
  px    <- mean(Slice(pxyz[, 1], begin))
  py    <- mean(Slice(pxyz[, 2], begin))
  pz    <- mean(Slice(pxyz[, 3], begin))
  stens <- mean(Slice(pxyz[, 2] - (pxyz[, 1] + pxyz[, 3]) * 0.5, begin))
  return (c(px, py, pz, stens))
}

observe.beg <- 100000
print(CalcPressure(LoadPressure("./box27.0x27.0"), observe.beg))
print(CalcPressure(LoadPressure("./box27.5x27.5"), observe.beg))
print(CalcPressure(LoadPressure("./box28.0x28.0"), observe.beg))
print(CalcPressure(LoadPressure("./box28.5x28.5"), observe.beg))
print(CalcPressure(LoadPressure("./box29.0x29.0"), observe.beg))
