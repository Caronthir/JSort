push!(LOAD_PATH, "$(homedir())/Development/")
using JSort

path = "$(homedir())/master/sortering/"*"28si.yaml"
parameters = Parameters(path)
#sortfile(parameters)
#savecoefficients_eΔe(parameters, coefficients...)
#@time events = loadlabr(parameters);
#@time makeeΔebin(events, parameters)
@time coeff_e, coeff_Δe = featurealign2d(parameters.readpath, plot=true)
@show coeff_e
p_e = ParticleCalibrator(coeff_e, var=:e)
p_Δe = ParticleCalibrator(coeff_e, var=:Δe)
write(joinpath(parameters.savepath, "coefficients_e.txt"), p_e)
write(joinpath(parameters.savepath, "coefficients_de.txt"), p_Δe)
#@time savecoefficients_eΔe(parameters, coefficients...)
