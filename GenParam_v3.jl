function GenParam(generators, setI, setI_hat, wind_LBset, wind_OSset, pvset, pv_BTMset, hydroset, nucset, coalset, oilset, gasset, σ0_pv, σ0_wind, υ0)
    filename_windcurve = string(pwd(),"/data/hourly_wind_gen_2018.csv")
    filename_pvcurve = string(pwd(),"/data/hourly_solar_gen_2018.csv")
    # pvcurve = GenCurveLoad(filename_pvcurve)[2:31:end,:]/2883.8  # pick one day a month over the year / nameplate capacity = 2883MW
    # windcurve = GenCurveLoad(filename_windcurve)[2:31:end,:]/1300 # pick one day a month over the year / nameplate capacity = 1300MW

    WTcf = 0.29
    PVcf = 0.14
    # https://www.nyiso.com/power-trends

    temp_pv = GenCurveLoad(filename_pvcurve)
    pvcurve_year = temp_pv / mean(temp_pv) * PVcf
    temp_wind = GenCurveLoad(filename_windcurve)
    windcurve_year = temp_wind / mean(temp_wind) * WTcf

    pvcurve = sum(pvcurve_year, dims=1)/365
    windcurve = sum(windcurve_year, dims=1)/365

    # DCRF = 6.3281e-4 # h=5years / 5%
    # DCRF = 3.5481e-4 * length(setT) / 24 # h=10years / 5%
    # DCRF = 2.1984e-4 # h=20years / 5%
    # DCRF = 3.05e-4 # h=10years / 2%
    # DCRF = 1.6755e-4 # h=20years / 2%

    Cinv = Dict()
    # Investment cost source:  https://www.eia.gov/electricity/generatorcosts/
    # PV & WT ref: https://www.eia.gov/outlooks/aeo/section_issue_renewables.php

    Cinv["Wind_LB"] = DCRF * 1260 * 1e3
    Cinv["Wind_OS"] = DCRF * 5446 * 1e3
    Cinv["PV"] = DCRF * 1307 * 1e3
    Cinv["PV_BTM"] = DCRF * 2019 * 1e3
    # Cinv["PV_BTM"] = DCRF * 1539 * 1e3
    # Cinv["PV_BTM"] = DCRF * 1307 * 1e3
    Cinv["Hydro"] = DCRF * 5312 * 1e3
    Cinv["Nuclear"] = DCRF * 1000 * 1e3
    Cinv["Coal"] = DCRF * 1000 * 1e3
    Cinv["Oil"] = DCRF * 1672 * 1e3
    # Cinv["Gas"] = DCRF * 1563 * 1e3
    Cinv["Gas"] = DCRF * 895 * 1e3

    Coper = Dict()
    Coper["Wind_LB"] = [0, 1.1, 0]
    Coper["Wind_OS"] = [0, 1.1, 0]
    Coper["PV"] = [0, 0.4, 0]
    Coper["PV_BTM"] = [0, 0.8, 0]
    # Coper["PV_BTM"] = [0, 0, 0]
    Coper["Hydro"] = [0, 4, 0]
    Coper["Nuclear"] = [0, 7, 1250]
    Coper["Coal"] = [0.0, 19.98, 1535.0]
    Coper["Oil"] = [0.0, 192.06, 5000.0]
    Coper["Gas"] = [0.0, 44.92, 612.84] # 20

    Cinvg = Vector{Float64}(undef, length(setI))
    Coperg = Vector{Vector{Float64}}(undef, length(setI))
    for g in setI_hat
        if g in wind_LBset
            Cinvg[g] = Cinv["Wind_LB"]
            Coperg[g] = Coper["Wind_LB"]
        elseif g in wind_OSset
            Cinvg[g] = Cinv["Wind_OS"]
            Coperg[g] = Coper["Wind_OS"]
        elseif g in pvset
            Cinvg[g] = Cinv["PV"]
            Coperg[g] = Coper["PV"]
        elseif g in pv_BTMset
            Cinvg[g] = Cinv["PV_BTM"]
            Coperg[g] = Coper["PV_BTM"]
        elseif g in hydroset
            Cinvg[g] = Cinv["Hydro"]
            Coperg[g] = Coper["Hydro"]
        elseif g in nucset
            Cinvg[g] = Cinv["Nuclear"]
            Coperg[g] = Coper["Nuclear"]
        elseif g in coalset
            Cinvg[g] = Cinv["Coal"]
            Coperg[g] = Coper["Coal"]
        elseif g in oilset
            Cinvg[g] = Cinv["Oil"]
            Coperg[g] = Coper["Oil"]
        elseif g in gasset
            Cinvg[g] = Cinv["Gas"]
            Coperg[g] = Coper["Gas"]
        end
    end


    for g in setdiff(setI, setI_hat)
        if g in wind_LBset
            Cinvg[g] = Cinv["Wind_LB"]
            Coperg[g] = generators[g].cost
        elseif g in wind_OSset
            Cinvg[g] = Cinv["Wind_OS"]
            Coperg[g] = generators[g].cost
        elseif g in pvset
            Cinvg[g] = Cinv["PV"]
            Coperg[g] = generators[g].cost
        elseif g in pv_BTMset
            Cinvg[g] = Cinv["PV_BTM"]
            Coperg[g] = generators[g].cost
        elseif g in hydroset
            Cinvg[g] = Cinv["Hydro"]
            Coperg[g] = generators[g].cost
        elseif g in nucset
            Cinvg[g] = Cinv["Nuclear"]
            Coperg[g] = generators[g].cost
        elseif g in coalset
            Cinvg[g] = Cinv["Coal"]
            Coperg[g] = generators[g].cost
        elseif g in oilset
            Cinvg[g] = Cinv["Oil"]
            Coperg[g] = generators[g].cost
        elseif g in gasset
            Cinvg[g] = Cinv["Gas"]
            Coperg[g] = generators[g].cost
        else
            Coperg[g] = generators[g].cost
        end
    end

    βFC = 1.0*ones(length(setI),setT[end])
    for g in union(pvset, pv_BTMset)
        for t in setT
            βFC[g,t] = pvcurve[1,t]
        end
    end
    for g in union(wind_LBset, wind_OSset)
        for t in setT
            βFC[g,t] = windcurve[1,t]
        end
    end

    σ = zeros(length(setI))
    for g in union(pvset, pv_BTMset)
        σ[g] = σ0_pv
    end
    for g in union(wind_LBset, wind_OSset)
        σ[g] = σ0_wind
    end

    υ = υ0*ones(length(setI), setT[end])
    # kappa = [1; 0.252; 0.75; 0.411; 0.44; 0.385]
    # PolicyB = 1e8  * ones(length(stateset)) # Policy Budget
    return Cinvg, Coperg, βFC, σ, υ
end
