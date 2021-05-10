using DelimitedFiles
function Demand24Load(filename_Load24)
    ## source from: http://www.pjm.com/markets-and-operations/ops-analysis/historical-load-data.aspx
    temp_data = readdlm(filename_Load24, ',',header=true)[1][:,2:end]
    raw_data = Array{Float64}(undef,size(temp_data))
    normalized_data = Array{Float64}(undef,size(temp_data))
    raw_data[:] = temp_data[:]
    max_load = [2586,1917,2744,538,1295,2358,2216,630,1388,9397,5214]
    N_bus = size(max_load,1)
    N_hours = size(temp_data,1)

    for k = 1:N_bus
        for m = 1:N_hours
        normalized_data[m,k] = raw_data[m,k] / max_load[k,1]
    end
end
    load24 = normalized_data
    return load24
end
