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
  p.tot <- mean(slice((pxyz[,1] + pxyz[,2] + pxyz[,3]) / 3.0, begin))
  return (c(px, py, pz, p.tot))
}

observe.beg <- 30000
print(calc.pressure(load.pressure("./exclude_myself/mol_dens20.0"), observe.beg))
print(calc.pressure(load.pressure("./exclude_myself/mol_dens30.0"), observe.beg))
print(calc.pressure(load.pressure("./exclude_myself/mol_dens35.0"), observe.beg))
print(calc.pressure(load.pressure("./exclude_myself/mol_dens40.0"), observe.beg))
print(calc.pressure(load.pressure("./exclude_myself/mol_dens45.0"), observe.beg))
print(calc.pressure(load.pressure("./exclude_myself/mol_dens50.0"), observe.beg))
print(calc.pressure(load.pressure("./exclude_myself/mol_dens60.0"), observe.beg))
