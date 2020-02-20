using .JSort
import PyCall

nndcurl(x) = "https://www.nndc.bnl.gov/nudat2/getdatasetClassic.jsp?nucleus=$(uppercase(x))&unc=nds"

function fetchlevels(isotope::AbstractString)
    #const url = nndcurl(isotope)
    #PyCall
    # c1 = GammaCascade(Gate(1779.03), [1778.969], [(1600, 1850)],
    #                   Dict(:align => Dict(:a1 => Dict(
    #                       :nregions=>2, :swidth=>60, :rwidth=>200, :plot=>true),
    #                                       :a2 => Dict(
    #                                           :nregions=>3, :swidth=>-10, :rwidth=>200, :plot=>true))))
    # c2 = GammaCascade(Gate(4617.86), [2838.29],  [(3000, 3400)],
    #                   Dict(:align => Dict(:a1 => Dict(
    #                       :nregions=>2, :swidth=>60, :rwidth=>200, :plot=>true),
    #                                       :a2 => Dict(
    #                                           :nregions=>3, :swidth=>-10, :rwidth=>200, :plot=>true))))
    # c3 = GammaCascade(Gate(4979.92), [3200.7],   [(3400, 4000)],
    #                   Dict(:align => Dict(:a1 => Dict(
    #                       :nregions=>2, :swidth=>60, :rwidth=>200, :plot=>true),
    #                                       :a2 => Dict(
    #                                           :nregions=>3, :swidth=>-10, :rwidth=>200, :plot=>true))))
    # c4 = GammaCascade(Gate(6276.20), [4496.92],  [(5400, 5700)],
    #                   Dict(:align => Dict(:a1 => Dict(
    #                       :nregions=>2, :swidth=>60, :rwidth=>200, :plot=>true),
    #                                       :a2 => Dict(
    #                                           :nregions=>3, :swidth=>-10, :rwidth=>200, :plot=>true))))
    # # WTF foregår her??
    # c5 = GammaCascade(Gate(6690.74, width=100), [4910.8],   [(6000, 6500)],
    #                   Dict(:align => Dict(:a1 => Dict(
    #                       :delta=>[2000, 2000], :nregions=>2, :swidth=>100, :rwidth=>300, :plot=>true),
    #                                       :a2 => Dict(
    #                       :delta=>[2000, 2000], :nregions=>3, :swidth=>-10, :rwidth=>200, :plot=>true),
    #                                       :bins2 => 2000)))
    # # Spør
    # c6 = GammaCascade(Gate(7799.01, width=70), [6018],   [(7500, 8000)],
    #                   Dict(:align => Dict(:a1 => Dict(
    #                       :delta=>[1000, 4000], :nregions=>2, :swidth=>450, :rwidth=>200, :plot=>true),
    #                                       :a2 => Dict(
    #                                           :nregions=>3, :swidth=>-10, :rwidth=>200, :plot=>true),
    #                                       ),
    #                        :bins2 => 2500))

    # c9 = GammaCascade(Gate(8904.4, width=200), [8901.8, 7123.8], [(8500, 9300), (7800, 8400)])
    # #extend!(c2, c1)
    #extend!(c3, [c2, c1])


    # Hvorfor forsvinner ikke grunntilstanden
    # Nydelig peak, men hva er de høyere peaksa?
    # De høyere forsvinner ved trangere gate
    c1 = GammaCascade(Gate2D(11033.623860362739, 11593.737177077102, 839.3651144773087, 1091.0419622932445),
                      [1778.969], [(1600, 1850)],
                      Dict(:align => Dict(
                          :a1 => Dict(
                              :nregions=>2, :swidth=>60, :rwidth=>200, :plot=>true),
                          :a2 => Dict(
                              :nregions=>3, :swidth=>-10, :rwidth=>200, :plot=>true))))
    # 1st of 16O? Gives two peaks. Probably not. 27Al?
    c2 = GammaCascade(Gate2D(11661.383712912168, 12108.019965763204, 854.5382927204712, 1011.440500564305),
                      [1778.969], [(1500, 1850)],
                      Dict(:align => Dict(
                          :a1 => Dict(
                              :nregions=>2, :swidth=>60, :rwidth=>200, :plot=>true),
                          :a2 => Dict(
                              :nregions=>3, :swidth=>-10, :rwidth=>200, :plot=>true))))
    c1
    #c1, c2, c3, c4, c5, c6, c9
end

