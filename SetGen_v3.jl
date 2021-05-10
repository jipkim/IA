function SetGen(generators, nodes, bus_gep)
   setN = 1:length(nodes)
   setI = 1:length(generators)

   wind_LBset = Int64[]
   wind_OSset = Int64[]
   pvset = Int64[]
   pv_BTMset = Int64[]
   hydroset = Int64[]
   nucset = Int64[]
   coalset = Int64[]
   oilset = Int64[]
   gasset = Int64[]

   for g in setI
      if generators[g].gtype == "Wind_LB"
         push!(wind_LBset, generators[g].gindex)
      elseif generators[g].gtype == "Wind_OS"
         push!(wind_OSset, generators[g].gindex)
      elseif generators[g].gtype == "PV"
         push!(pvset, generators[g].gindex)
      elseif generators[g].gtype == "PV_BTM"
         push!(pv_BTMset, generators[g].gindex)
      elseif generators[g].gtype == "Hydro"
         push!(hydroset, generators[g].gindex)
      elseif generators[g].gtype == "Nuclear"
         push!(nucset, generators[g].gindex)
      elseif generators[g].gtype == "Coal"
         push!(coalset, generators[g].gindex)
      elseif generators[g].gtype == "Oil"
         push!(oilset, generators[g].gindex)
      elseif generators[g].gtype == "Gas"
         push!(gasset, generators[g].gindex)
      end
   end

   # Avoid MIP (Not implementing UC problem)
   for g in setI
      generators[g].Pmin = 0
   end

   setI_hat = setI[end]+1 : setI[end]+5
   push!(wind_LBset, setI[end]+1)
   push!(wind_OSset, setI[end]+2)
   push!(pvset, setI[end]+3)
   push!(pv_BTMset, setI[end]+4)
   push!(gasset, setI[end]+5)

   B_g = []
   for g in setI
      push!(B_g, generators[g].location)
   end

   for jj in 1:5
      push!(B_g, bus_gep)
   end

   B_gn = Array{Array{Int64},1}(undef,length(nodes))
   for ii in 1:length(nodes)
      B_gn[ii] = Int64[]
   end
   for g in union(setI,setI_hat)
      push!(B_gn[B_g[g]], g)
   end

   setI = union(setI, setI_hat)

   I_R = intersect(setI, union(wind_LBset,wind_OSset,pvset,pv_BTMset,hydroset))
   I_R_hat = intersect(setI_hat, union(wind_LBset,wind_OSset,pvset,pv_BTMset,hydroset))
   I_C = setdiff(setI, union(wind_LBset,wind_OSset,pvset,pv_BTMset,hydroset))
   I_C_hat = setdiff(setI_hat,union(wind_LBset,wind_OSset,pvset,pv_BTMset,hydroset))

   return wind_LBset, wind_OSset, pvset, pv_BTMset, hydroset, nucset, coalset, oilset, gasset, B_g, B_gn, setI, setI_hat, I_R, I_R_hat, I_C, I_C_hat
end
