push!(LOAD_PATH, "$(homedir())/Development/")
using JSort

path = "$(homedir())/master/sortering/"*"zn70.yaml"
@time parameters = Parameters(path)
@time sortfile(parameters)
