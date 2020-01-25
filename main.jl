push!(LOAD_PATH, "$(homedir())/Development/")
using JSort

path = "$(homedir())/master/sortering/"
parameters = Parameters(path*"si28.yaml")
coefficients = featurealign2d("/home/erdos/master/sortering/sirius/")
savecoefficients_eΔe(parameters, coefficients...)
println(coefficients)
