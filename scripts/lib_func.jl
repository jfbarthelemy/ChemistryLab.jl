using SymPy

x, α = symbols("x α", positive=true)

diff(x^α, x)
diff(x^(0), x)
diff(x^(-1), x)
diff(x^(-2), x)
diff(x^(-3), x)
diff(log(x), x)
diff(log(x)/x, x)

integrate(x^α, x)
integrate(x^(-1), x)
integrate(x^(-2), x)
integrate(x^(-3), x)
integrate(log(x), x)
integrate(log(x)/x, x)

integrate(integrate(x^α, x), x)
integrate(integrate(x^(-1), x), x)
integrate(integrate(x^(-2), x), x)
integrate(integrate(x^(-3), x), x)
integrate(integrate(log(x), x), x)
integrate(integrate(log(x)/x, x), x)
