#dat <- data.matrix(read.table("./angle_hist.txt"))
#dat <- data.matrix(read.table("./h2e.txt"))
dat <- data.matrix(read.table("./bondht_hist.txt"))
mean(dat)
max(dat)

x <- hist(dat, breaks = 100)


plot(x$mids, x$density)
plot(x$mids, x$density / sin(x$mids))
