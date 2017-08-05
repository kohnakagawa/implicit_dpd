kappa = 100.0
rho_co = 40.0
Reo = 3.5
v = -2.0 * (kappa + 3.0) / rho_co
w = (3.0 / 2.0) * (kappa + 2.0) / (rho_co ** 2)
f(x) = (x + v * 0.5 * x**2 + (2.0 * w / 3.0) * x**3) / Reo**3
plot "rho_vs_p.dat"
replot f(x)