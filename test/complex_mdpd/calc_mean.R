LoadPressure <- function(root.dir) {
  fname <- file.path(root.dir, "pressure.txt")
  pxyz  <- data.matrix(read.table(fname))
  return (pxyz)
}

Slice <- function(vec, begin) {
  return (vec[seq(begin, length(vec))])
}

CalcPressure <- function(pxyz, begin, Ly) {
  px    <- mean(Slice(pxyz[, 1], begin))
  py    <- mean(Slice(pxyz[, 2], begin))
  pz    <- mean(Slice(pxyz[, 3], begin))
  stens <- mean(Slice(pxyz[, 2] - (pxyz[, 1] + pxyz[, 3]) * 0.5, begin))
  print(paste("Ly=", Ly, ": p = c(", px, py, pz, ")", sep=" ") )
  return (stens * Ly)
}

stens.mean <- c()
observe.beg <- 100000
box.leng <- c(27.25, 27.5, 27.75, 28.0, 28.25, 28.5, 28.75, 29.0)
root.dir <- "./exclude_myself"
for (blen in box.leng) {
  fname <- file.path(root.dir, paste("box", blen, "x", blen, sep=""))
  stens.mean <- c(stens.mean, CalcPressure(LoadPressure(fname), observe.beg, blen))
}

print(stens.mean)

area <- box.leng * box.leng
dat <- data.frame(area, stens.mean)
print(dat)
lm(stens.mean ~ area, data = dat)
plot(dat)
