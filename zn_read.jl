push!(LOAD_PATH, "$(homedir())/Development/")
using JSort

path = "$(homedir())/master/sortering/"*"zn70.yaml"
path = "$(homedir())/master/sortering/28si.yaml"
@time parameters = Parameters(path)
@time reader = Reader(parameters)
@time read(reader)
