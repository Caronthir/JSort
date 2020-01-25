using .JSort

function sortfile(parameters::Parameters)
    labr = Vector{LaBrEvent}
    for file in parameters.datafiles
        println("Processing file $file")
        @time labr = processfile(file, parameters)
    end
    println("Saving.")
    @show length(labr)
    @time save(labr, parameters.savepath)
end
