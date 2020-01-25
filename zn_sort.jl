push!(LOAD_PATH, "$(homedir())/Development/")
using JSort

path = "$(homedir())/master/sortering/"*"zn70.yaml"
parameters = Parameters(path)
#sortfile(parameters)
#savecoefficients_eΔe(parameters, coefficients...)
#@time events = loadlabr(parameters);
#@time makeeΔebin(events, parameters)
@time coefficients = featurealign2d("/home/erdos/master/sortering/zinc", plot=true)
@time savecoefficients_eΔe(parameters, coefficients...)
