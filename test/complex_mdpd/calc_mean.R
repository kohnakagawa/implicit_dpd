load.pressure <- function(root.dir) {
  fname <- file.path(root.dir, "pressure.txt")
  pxyz <- data.matrix(read.table(fname))
  return (pxyz)
}

slice <- function(vec, begin) {
  return (vec[seq(begin, length(vec))])
}

calc.pressure <- function(pxyz, begin) {
  px <- mean(slice(pxyz[,1], begin))
  py <- mean(slice(pxyz[,2], begin))
  pz <- mean(slice(pxyz[,3], begin))
  stens <- mean(slice(pxyz[,2] - (pxyz[,1] + pxyz[,3]) * 0.5, begin))
  return (c(px, py, pz, stens))
}

observe.beg <- 100000
print(calc.pressure(load.pressure("./box27.0x27.0"), observe.beg))
print(calc.pressure(load.pressure("./box27.5x27.5"), observe.beg))
print(calc.pressure(load.pressure("./box28.0x28.0"), observe.beg))
print(calc.pressure(load.pressure("./box28.5x28.5"), observe.beg))
print(calc.pressure(load.pressure("./box29.0x29.0"), observe.beg))
