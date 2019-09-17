using .JSort
using Distributed
using SharedArrays

function sortfile(parameterpath::T) where T<:AbstractString
    parameters = Parameters(parameterpath)
    if parameters.cores > 1
        println("Using $(parameters.cores) cores.")
        addprocs(parameters.cores)
        # @everywhere using SharedArrays
    end

    println("Processing files.")
    labr = Vector{LaBrEvent}
    for file in parameters.datafiles
        @time labr = processfile(file, parameters)
    end
    println("Saving.")
    #@show sum(matrices[:eΔe].matrix)
    @show length(labr)
    #save(matrices, "sirius")
    save(labr, parameters.savepath)
end
