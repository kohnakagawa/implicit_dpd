#!/usr/bin/gnuplot

set print 'para0.dat'
list = system('ls sample*/meanspect.txt')
ranges = "0.75 1.0 1.5 2.0 2.5 3.0"
do for [ i = 1 : words(list)] {
    do for [u = 1 : words(ranges)] {
        EXPR = "f(x) = a*x + b*x**2"
        eval(EXPR)
        FIT = sprintf( "fit [:%s] f(x) '%s' u ($1*$1):(1/$3) via a,b", word(ranges, u), word(list, i))
        eval(FIT)
        print word(ranges, u), " ", b
    }
    print "\n"
}