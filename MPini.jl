using Pkg, DataFrames, Distributions, DelimitedFiles, CSV
using JuMP, Gurobi, Ipopt
using KNITRO
function MPini(testsystemT, testsystemD, busT_gep, busD_gep, κ, setT, η, σ0_pv, σ0_wind, υ0, Solver, bounds::bounds, D_ratio, ScaleObj, SD_tol, M, N, ED, γCT, ST_res, WMP, CTI)   

   include("NetworkDataType.jl")
   include("NetworkLoad_v3.jl")   
   include("Demand24Load_NYISO.jl")
   filename_Load24 = string(pwd(), "/data/load24data-NYISO.csv")
   temp_loadcurve = Demand24Load(filename_Load24)
   loadcurveT = D_scale * temp_loadcurve'   
   loadcurveD = D_scale * [0.5915717782270937, 0.5603916143450037, 0.5383633074385442, 0.524103437267213, 0.5188890071299351, 0.5265510269234862, 0.5592210279876556, 0.6044482281579228, 0.6518037671597319, 0.6848994359902096, 0.7139512610407577, 0.7381079067787591, 0.7567308715547515, 0.7688624028945408, 0.7789720123443652, 0.7858891135468766, 0.7866340321379164, 0.7777003298925189, 0.7514100244758966, 0.7289560498031287, 0.7203362775353837, 0.7102266680855592, 0.6857507715228265, 0.6504203469192296]

   busesT, linesT, generatorsT, datamatT, setL_T, setN_T, DTp, DTq = NetworkLoad(testsystemT)
   busesD, linesD, generatorsD, datamatD, setL_D, setN_D, DDp, DDq = NetworkLoad(testsystemD)
   Γpf = map((x) -> if isnan(x) 0 else x end, DDq./DDp)

   
   include("SetGen_v3.jl")
   wind_LBsetT, wind_OSsetT, pvsetT, pv_BTMsetT, hydrosetT, nucsetT, coalsetT, oilsetT, gassetT, B_gT, B_gnT, setI_TS, setI_TS_hat, setI_RT, setI_RT_hat, setI_CT, setI_CT_hat = SetGen(generatorsT, busesT, busT_gep)
   wind_LBsetD, wind_OSsetD, pvsetD, pv_BTMsetD, hydrosetD, nucsetD, coalsetD, oilsetD, gassetD, B_gD, B_gnD, setI_DS, setI_DS_hat, setI_RD, setI_RD_hat, setI_CD, setI_CD_hat = SetGen(generatorsD, busesD, busD_gep)

   include("GenCurveLoad.jl")
   include("GenParam_v3.jl")
   CinvgT, CopergT, βFC_T, σT, υT = GenParam(generatorsT, setI_TS, setI_TS_hat, wind_LBsetT, wind_OSsetT, pvsetT, pv_BTMsetT, hydrosetT, nucsetT, coalsetT, oilsetT, gassetT, σ0_pv, σ0_wind, υ0)
   CinvgD, CopergD, βFC_D, σD, υD = GenParam(generatorsD, setI_DS, setI_DS_hat, wind_LBsetD, wind_OSsetD, pvsetD, pv_BTMsetD, hydrosetD, nucsetD, coalsetD, oilsetD, gassetD, σ0_pv, σ0_wind, υ0)

   wind_LBsetD = map((x) -> x + 1000, wind_LBsetD)
   wind_OSsetD = map((x) -> x + 1000, wind_OSsetD)
   pvsetD = map((x) -> x + 1000, pvsetD)
   pv_BTMsetD = map((x) -> x + 1000, pv_BTMsetD)
   hydrosetD = map((x) -> x + 1000, hydrosetD)
   nucsetD = map((x) -> x + 1000, nucsetD)
   coalsetD = map((x) -> x + 1000, coalsetD)
   oilsetD = map((x) -> x + 1000, oilsetD)
   gassetD = map((x) -> x + 1000, gassetD)
   setI_DS = map((x) -> x + 1000, setI_DS)
   setI_DS_hat = map((x) -> x + 1000, setI_DS_hat)
   setI_RD = map((x) -> x + 1000, setI_RD)
   setI_RD_hat = map((x) -> x + 1000, setI_RD_hat)
   setI_CD = map((x) -> x + 1000, setI_CD)
   setI_CD_hat = map((x) -> x + 1000, setI_CD_hat)

   Pmax = merge(Dict(i => generatorsT[i].Pmax for i = 1:length(generatorsT)),
               Dict(i+1000 => generatorsD[i].Pmax for i = 1:length(generatorsD)))
   Qmax = merge(Dict(i => generatorsT[i].Qmax for i = 1:length(generatorsT)),
               Dict(i+1000 => generatorsD[i].Qmax for i = 1:length(generatorsD)))

   B_gT = Dict(i => B_gT[i] for i = 1:length(B_gT))
   B_gD = Dict(i + 1000 => B_gD[i] for i = 1:length(B_gD))
   B_g = merge(B_gT, B_gD)
   B_gnD = merge(Dict(i => B_gnD[i] + 1000 * ones(Int64, length(B_gnD[i]),1) for i = 1:length(B_gnD) if ~isempty(B_gnD[i])), 
      Dict(i => B_gnD[i] for i = 1:length(B_gnD) if isempty(B_gnD[i])))

   setI_R = union(setI_RT, setI_RD)
   setI_R_hat = union(setI_RT_hat, setI_RD_hat)
   setI_C = union(setI_CT, setI_CD)
   setI_C_hat = union(setI_CT_hat, setI_CD_hat)
   setI_T = union(setI_TS, setI_DS)
   setI_T_hat = union(setI_TS_hat, setI_DS_hat)
   setI_D = setI_DS
   setI_D_hat = setI_DS_hat
   setI_DER = union(setI_DS, B_gnT[busT_gep])
   setI_DER_hat = union(setI_DS_hat, intersect(setI_T_hat, B_gnT[busT_gep]))

   CCT = 50
   H = Dict(i =>  if i in union(gassetT, gassetD) 0.419
         elseif i in union(coalsetT, coalsetD) 1.002
         elseif i in union(oilsetT, oilsetD) 0.9606 
         else 0 end
      for i in setI_C)

   Cinvg = merge(Dict(i => CinvgT[i] for i = 1:length(CinvgT)), Dict(i + 1000 => CinvgD[i] for i = 1:length(CinvgD)))
   Coperg = merge(Dict(i => CopergT[i] for i = 1:length(CopergT)), Dict(i + 1000 => CopergD[i] for i = 1:length(CopergD)))
   βFC = merge(Dict([i,t] => βFC_T[i,t] for i in 1:size(βFC_T,1), t in setT), 
               Dict([i + 1000, t] => βFC_D[i,t] for i in 1:size(βFC_D,1), t in setT))
   σ = merge(Dict(i => σT[i] for i = 1:length(σT)), Dict(i + 1000 => σD[i] for i = 1:length(σD)))
   υ = merge(Dict([i,t] => υT[i,t] for i = 1:size(υT,1), t in setT), 
            Dict([i + 1000,t] => υD[i,t] for i = 1:size(υD,1), t in setT))

   m = Model(KNITRO.Optimizer)
   set_optimizer_attribute(m, "ms_enable", 1)
   set_optimizer_attribute(m, "ms_maxsolves", 300)
   set_optimizer_attribute(m, "ms_terminate", 0)
   set_optimizer_attribute(m, "par_blasnumthreads", 10)
   set_optimizer_attribute(m, "maxit", 1000000)                     

   @variables(m, begin
      bounds.πC_LB <= πC[n in setN_D, t in setT] <= bounds.πC_UB
      bounds.πDER_LB <= πDER[i in setI_D, t in setT] <= bounds.πDER_UB
      bounds.s_LB <= s[n in setN_D, t in setT] <= bounds.s_UB
      bounds.τe_LB <= τe <= bounds.τe_UB
      bounds.τc_LB <= τc <= bounds.τc_UB
      Obj_REG   
      0 <= gTpo[i in setI_T, t in setT]
      0 <= gDpo[i in setI_D, t in setT]
      0 <= gDqo[i in setI_D, t in setT]
      0 <= gpmax[i in setI_T_hat]
      Obj_DER
   end)
   @constraint(m, gpmax_fix1[i in intersect(pv_BTMsetT, setI_T_hat)], gpmax[i] == 0)
   @constraint(m, gpmax_fix2[i in intersect(wind_LBsetT, setI_T_hat)], gpmax[i] == 0)
   @constraint(m, gpmax_fix3[i in intersect(wind_OSsetD, setI_T_hat)], gpmax[i] == 0)
   @constraint(m, gpmax_fix4[i in intersect(wind_LBsetD, setI_T_hat)], gpmax[i] == 0)
   @constraint(m, gpmax_fix5[i in intersect(pvsetD, setI_T_hat)], gpmax[i] == 0)
   @constraint(m, gpmax_fix6[i in intersect(gassetD, setI_T_hat)], gpmax[i] == 0)
   @variables(m, begin   
      0 <= gDp[i in setI_D, t in setT]
      0 <= gDq[i in setI_D, t in setT]
      0 <= aD[l in setL_D, t in setT]
      0 <= vD[n in setN_D, t in setT]
      fDp[l in setL_D, t in setT]
      fDq[l in setL_D, t in setT]
      Obj_D_prim[t in setT]
      Obj_D_prim_act[t in setT]
      OCgenD[i in setI_D, t in setT]
      Obj_C[t in setT]
      fTp[l in setL_T, t in setT]
      0 <= gTp[i in setI_T, t in setT]
      θT[n in setN_T, t in setT]
      Obj_W_prim[t in setT]      
      OCgenT[i in setI_T, t in setT]
      ξT[l in setL_T, t in setT]
      λT[setn in setN_T, t in setT]
      0 <= γTlb[i in setI_T, t in setT]
      0 <= γTub[i in setI_T, t in setT]
      0 <= δTlb[l in setL_T, t in setT]
      0 <= δTub[l in setL_T, t in setT]
      Obj_W_dual[t in setT]    
      λD[n in setN_D, t in setT]
      μD[n in setN_D, t in setT]
      βD[l in setL_D, t in setT]            
      0 <= ζDa[l in setL_D, t in setT]
      ζDs[l in setL_D, t in setT]
      ζDp[l in setL_D, t in setT]
      ζDq[l in setL_D, t in setT]
      0 <= ηDsub[l in setL_D, t in setT]
      ηDpub[l in setL_D, t in setT]
      ηDqub[l in setL_D, t in setT]
      0 <= ηDslb[l in setL_D, t in setT]
      ηDplb[l in setL_D, t in setT]
      ηDqlb[l in setL_D, t in setT]
      0 <= γDlb[i in setI_D, t in setT]
      0 <= γDub[i in setI_D, t in setT]
      0 <= δDlb[i in setI_D, t in setT]
      0 <= δDub[i in setI_D, t in setT]
      0 <= ϑDlb[n in setN_D, t in setT]
      0 <= ϑDub[n in setN_D, t in setT]
      Obj_D_dual[t in setT]
   end)
   d = DDp * loadcurveD'
   @variable(m, MPbound)
   @constraint(m, MPbound == Obj_REG)
   @objective(m, Max, Obj_REG)
   @constraint(m, Obj_REG_def, Obj_REG == 
         Obj_DER         
         + sum(Obj_D_prim_act[t] + Obj_C[t] + Obj_W_prim[t] for t in setT)
         - sum(ED * H[i] * gTp[i,t] for i in setI_C, t in setT)
         - sum(ED * H[i] * gDp[i,t] for i in intersect(setI_C, setI_D), t in setT)
         + (1-γCT) * (
               sum(CCT * H[i] * gTp[i,t] for i in setdiff(setI_C, setI_D), t in setT)                  
            ) 
      )

   @constraint(m, Subsidy_budget, 
      sum(
            s[1,t] * (DTp[busT_gep] * loadcurveT[busT_gep,t]) 
            + sum(s[n,t] * d[n,t] for n in setN_D)
            + sum(ST_res[n,t] * DTp[n] * loadcurveT[n,t] for n in setdiff(setN_T, [busT_gep]))
            + γCT * (
               sum(CCT * H[i] * gTp[i,t] for i in setdiff(setI_C, setI_D))                   
            ) 
         for t in setT)
      >= 
      sum(
            τe * sum(gTp[i,t] for i in setdiff(setI_R, setI_D))            
         for t in setT)
      + τc * DCRF * sum(gpmax[i] for i in setdiff(setI_R_hat, setI_D))
   )
  
   F = 9.13
   hM = 0.419
   @variable(m, πDER_single)
   @variable(m, πC_single)
   @variable(m, s_single)
   @constraint(m, s_single_c[n in setN_D, t in setT], s[n,t] == s_single)

   if policy_switch == 1
      @constraint(m, piC_single[n in setN_D, t in setT], πC[n,t] == πC_single)
      @constraint(m, piDER_single[i in setI_D, t in setT], πDER[i,t] == πDER_single)
      @constraint(m, piEqual, πDER_single == πC_single)
   elseif policy_switch == 2
      @constraint(m, piC_single[n in setN_D, t in setT], πC[n,t] == πC_single)
      @constraint(m, piDER_VS_r[i in intersect(setI_D, setI_R), t in setT], πDER[i,t] == λT[busT_gep,t] .+ hM*ED .+ F)
      @constraint(m, piDER_VS_c[i in intersect(setI_D, setI_C), t in setT], πDER[i,t] == λT[busT_gep,t] .+ F)           
   elseif policy_switch == 4
      @constraint(m, piC_DLMP[n in setN_D, t in setT], πC[n,t] == λD[n,t])
      @constraint(m, piDER_DLMP[i in setI_D, t in setT], πDER[i,t] == λD[B_g[i],t])
   end

   @constraint(m, RPS, 
         sum(gTp[i,t] for i in setI_R, t in setT) 
         + sum(gDp[i,t] for i in intersect(setI_R, setI_D), t in setT)         
         >=
         κ * (
               sum(DTp[n] * loadcurveT[n,t] for n in setN_T, t in setT)
               + sum(d[n,t] for n in setN_D, t in setT) 
            )
      )
   @constraint(m, CR_DER, Obj_DER >= 0)
   if policy_switch != 4 @constraint(m, CR_Utility, sum(Obj_D_prim_act[t] for t in setT) >= 0) end
   if D_new == 1
      @constraint(m, Obj_C_def[t in setT], Obj_C[t] == sum(
         ifelse(d0[n,t]==0, 0, 
               M_new * d[n,t] 
               - 0.5 * N_new[n,t] * d[n,t]^2 
               - (πC[n,t] + s[n,t]) * d[n,t]
            )         
      for n in setN_D)
         + M_new * (DTp[busT_gep] * loadcurveT[busT_gep,t])  
         - 0.5 * N_new_DTp[t] * (DTp[busT_gep] * loadcurveT[busT_gep,t])^2
         - (πC[1,t] + s[1,t]) * (DTp[busT_gep] * loadcurveT[busT_gep,t]) 
         )      
   else
      @constraint(m, Obj_C_def[t in setT], Obj_C[t] == sum(
         M[n,t] * d[n,t] 
         - 0.5 * N * d[n,t]^2 
         - (πC[n,t] + s[n,t]) * d[n,t]
      for n in setN_D)
         + M[1,t] * (DTp[busT_gep] * loadcurveT[busT_gep,t])  
         - 0.5 * N * (DTp[busT_gep] * loadcurveT[busT_gep,t])^2
         - (πC[1,t] + s[1,t]) * (DTp[busT_gep] * loadcurveT[busT_gep,t]) 
         )
   end

   

   if IA_switch == 1
      @constraint(m, Obj_DER_def, Obj_DER == sum(
                  sum((λT[busT_gep,t] - Coperg[i][2]) * gTp[i,t] for i in setI_DER)                
                  + sum((πDER[i,t] - Coperg[i][2]) * gDp[i,t] for i in intersect(setI_DER, setI_D))
                  + sum(τe * gTp[i,t] for i in setdiff(intersect(setI_DER, setI_R), setI_D))                  
                  - sum(CCT * H[i] * gTp[i,t] for i in setdiff(intersect(setI_DER, setI_C), setI_D))
         for t in setT)
         + sum(τc * DCRF * gpmax[i] for i in intersect(setI_DER_hat, setI_R))
         - sum(Cinvg[i] * gpmax[i] for i in setI_DER_hat)
         )
      
   else    
      K_VaR = quantile.(Normal(), 1-ϵ)
      K_CVaR = pdf(Normal(), quantile.(Normal(), 1-ϵ)) / (ϵ)    
      @variables(m, begin
               zAux[i in intersect(setI_DER, setI_D), t in setT]
               μDER[i in intersect(setI_DER, setI_D), t in setT]
               σDER[i in intersect(setI_DER, setI_D), t in setT]
         end)
      @constraint(m, Obj_DER_def, Obj_DER == sum(
                  sum((λT[busT_gep,t] - Coperg[i][2]) * gTp[i,t] for i in setI_DER)                
                  + sum((- Coperg[i][2]) * gDp[i,t] for i in intersect(setI_DER, setI_D))
                  + sum(zAux[i,t] for i in intersect(setI_DER, setI_D))
                  + sum(τe * gTp[i,t] for i in setdiff(intersect(setI_DER, setI_R), setI_D))                  
                  - sum(CCT * H[i] * gTp[i,t] for i in setdiff(intersect(setI_DER, setI_C), setI_D))
         for t in setT)
         + sum(τc * DCRF * gpmax[i] for i in intersect(setI_DER_hat, setI_R))
         - sum(Cinvg[i] * gpmax[i] for i in setI_DER_hat)
         )
      @constraint(m, zAux_reform[i in intersect(setI_DER, setI_D), t in setT], (μDER[i,t] - K_CVaR * σDER[i,t]) * gDp[i,t] >= zAux[i,t])    
      if policy_switch == 1
         @constraint(m, μDER_fix[i in intersect(setI_DER, setI_D), t in setT], μDER[i,t] == 1*πDER_single)
         @constraint(m, σDER_fix[i in intersect(setI_DER, setI_D), t in setT], σDER[i,t] == 0.15*μDER[i,t])             
      elseif policy_switch == 2
         @constraint(m, μDER_fix_r[i in intersect(setI_D, setI_R), t in setT], μDER[i,t] == 1*λT[busT_gep,t] .+ hM*ED .+ 1F)
         @constraint(m, μDER_fix_c[i in intersect(setI_D, setI_C), t in setT], μDER[i,t] == 1*λT[busT_gep,t] .+ 1F)
         @constraint(m, σDER_fix_r[i in intersect(setI_D, setI_R), t in setT], σDER[i,t] == 0.15*λT[busT_gep,t] .+ 0hM*ED .+ 0.15F)
         @constraint(m, σDER_fix_c[i in intersect(setI_D, setI_C), t in setT], σDER[i,t] == 0.15*λT[busT_gep,t] .+ 0.15F)
      elseif policy_switch == 4
         @constraint(m, μDER_fix[i in intersect(setI_DER, setI_D), t in setT], μDER[i,t] == 1*λD[B_g[i],t])
         @constraint(m, σDER_fix[i in intersect(setI_DER, setI_D), t in setT], σDER[i,t] == 0.15*μDER[i,t])        
      end
   end


   @constraint(m, WMP_bb[i in setdiff(setI_DER, B_gnT[busT_gep]), t in setT], gTpo[i,t] == 0)
   @constraint(m, gpo_UB_TMonly[i in intersect(setdiff(setI_DER, setI_DER_hat), B_gnT[busT_gep]), t in setT], gTpo[i,t] == βFC[[i,t]] * Pmax[i] + υ[[i,t]])
   @constraint(m, gpo_UB_both[i in setdiff(setdiff(setI_DER, setI_DER_hat), B_gnT[busT_gep]), t in setT], gTpo[i,t] + gDpo[i,t] == βFC[[i,t]] * Pmax[i] + υ[[i,t]])
   @constraint(m, gpo_UB_hat_TMonly[i in intersect(setI_DER_hat, B_gnT[busT_gep]), t in setT], gTpo[i,t] == βFC[[i,t]] * gpmax[i] + υ[[i,t]])
   @constraint(m, gpo_UB_hat_both[i in setdiff(setI_DER_hat, B_gnT[busT_gep]), t in setT], gTpo[i,t] + gDpo[i,t] == βFC[[i,t]] * gpmax[i] + υ[[i,t]])
   @constraint(m, gDqo_fix[i in setdiff(setI_D, setI_D_hat), t in setT], gDqo[i,t] == βFC[[i,t]] * Qmax[i] + υ[[i,t]])
   @constraint(m, gDqo_fix_hat[i in setI_D_hat, t in setT], gDqo[i,t] == βFC[[i,t]] * gpmax[i] + υ[[i,t]])
   @constraint(m, gTpo_fix[i in setdiff(setI_T, setI_DER), t in setT], gTpo[i,t] == βFC[[i,t]] * Pmax[i] + υ[[i,t]])   
   @variable(m, fDp0[t in setT])
   @constraint(m, fDp0_def[t in setT], fDp0[t] == D_ratio * sum(d[n,t] for n in setN_D))

   if policy_switch == 4
    @constraint(m, O_Dprim[t in setT], Obj_D_prim[t] == 
        - sum( (Coperg[i][2]) * gDp[i,t] for i in intersect(setI_D, setI_R))
        - sum( (Coperg[i][2]) * gDp[i,t] for i in intersect(setI_D, setI_C))
        )
    @constraint(m, O_Dprim_act_def[t in setT], Obj_D_prim_act[t] == 
            Obj_D_prim[t]            
        )
   else
      @constraint(m, O_Dprim[t in setT], Obj_D_prim[t] == - sum(πDER[i,t] * gDp[i,t] for i in setI_D))
      @constraint(m, O_Dprim_act_def[t in setT], Obj_D_prim_act[t] == 
         + sum(πC[n,t] * d[n,t] for n in setN_D)
         - sum(πDER[i,t] * gDp[i,t] for i in setI_D)
         - λT[busT_gep,t] * (fDp0[t])         
         )
   end      
   @constraint(m, BetweenNodes[l in setL_D, t in setT], 
         1e6(
            vD[linesD[l].fbus,t] 
            - 2*(linesD[l].r*fDp[l,t] + linesD[l].x*fDq[l,t]) 
            + aD[l,t]*(linesD[l].r^2 + linesD[l].x^2)
         )
         == 1e6vD[linesD[l].tbus,t])
   
   @constraint(m, RootV[t in setT], vD[1,t] == 1)

   @constraint(m, PBalance_root[n in [1], t in setT], 
      sum(fDp[l,t] for l in busesD[n].outline)      
      - sum(- linesD[l].r*aD[l,t] for l in busesD[n].inline)  
      - (fDp0[t] - sum(gTp[i,t] for i in setI_D))           
      - sum(gTp[i,t] + gDp[i,t] for i in B_gnD[n])
      + d[n,t] + vD[n,t]*busesD[n].Gs == 0)

   @constraint(m, PBalance[n in setdiff(setN_D,[1]), t in setT], 
      sum(fDp[l,t] for l in busesD[n].outline)  
      - sum(fDp[l,t] - linesD[l].r*aD[l,t] for l in busesD[n].inline) 
      - sum(gTp[i,t] + gDp[i,t] for i in B_gnD[n])
      + d[n,t] + vD[n,t]*busesD[n].Gs == 0)

   @constraint(m, QBalance[n in setN_D, t in setT], 
      sum(fDq[l,t] for l in busesD[n].outline)  
      - sum(fDq[l,t] - linesD[l].x*aD[l,t] for l in busesD[n].inline) 
      - sum(gDq[i,t] for i in B_gnD[n]) 
      + Γpf[n] * d[n,t] - vD[n,t]*busesD[n].Bs == 0)


   @constraint(m, SOCP[l in setL_D, t in setT], (fDp[l,t])^2 + (fDq[l,t])^2 <= vD[linesD[l].fbus,t] * aD[l,t])
   @constraint(m, LineCapFW[l in setL_D, t in setT], fDp[l,t]^2 + fDq[l,t]^2 <= linesD[l].u^2)
   @constraint(m, LineCapBW[l in setL_D, t in setT], (fDp[l,t] - aD[l,t] * linesD[l].r)^2 + (fDq[l,t] - aD[l,t] * linesD[l].x)^2 <= linesD[l].u^2)
   
   @constraint(m, PGmax[i in setI_D, t in setT], gDp[i,t] <= gDpo[i,t])
   @constraint(m, QGmax[i in setI_D, t in setT], gDq[i,t] <= gDqo[i,t])
   @constraint(m, Vmax[n in setN_D, t in setT], vD[n,t] <= busesD[n].Vmax^2)
   @constraint(m, Vmin[n in setN_D, t in setT], vD[n,t] >= busesD[n].Vmin^2)
   @constraint(m, Obj_D_dual_def[t in setT], Obj_D_dual[t] ==         
         sum(- (d[n,t] - (fDp0[t] - sum(gTp[i,t] for i in setI_D)))* λD[n,t] - Γpf[n] * d[n,t] * μD[n,t] for n in [1])
         + sum(- d[n,t] * λD[n,t] - Γpf[n] * d[n,t] * μD[n,t] for n in setdiff(setN_D,[1]))
         + sum(- busesD[n].Vmin^2 * ϑDlb[n,t] + busesD[n].Vmax^2 * ϑDub[n,t] for n in setdiff(setN_D,1))
         + sum(- ϑDlb[n,t] + ϑDub[n,t] for n in [1])
         + sum(linesD[l].u * (ηDslb[l,t] + ηDsub[l,t]) for l in setL_D)
         + sum(gDpo[i,t] * γDub[i,t] + gDqo[i,t] * δDub[i,t] for i in setI_D)
         - sum(- gTp[i,t] * λD[B_gD[i],t] for i in setI_D)       
      )
   
   @constraint(m, L2dual_SOC_etaub[l in setL_D, t in setT], ηDpub[l,t]^2 + ηDqub[l,t]^2 <= ηDsub[l,t]^2)
   @constraint(m, L2dual_SOC_etalb[l in setL_D, t in setT], ηDplb[l,t]^2 + ηDqlb[l,t]^2 <= ηDslb[l,t]^2)
   @constraint(m, L2dual_SOC_zeta[l in setL_D, t in setT], ζDs[l,t]^2 + ζDp[l,t]^2 + ζDq[l,t]^2 <= ζDa[l,t]^2)

   if policy_switch == 4
      @constraint(m, L2dual_gDp_DLMP_R[i in intersect(setI_D, setI_R), t in setT], - (Coperg[i][2] ) + λD[B_gD[i],t] + γDlb[i,t] - γDub[i,t] == 0)
      @constraint(m, L2dual_gDp_DLMP_C[i in intersect(setI_D, setI_C), t in setT], - (Coperg[i][2] ) + λD[B_gD[i],t] + γDlb[i,t] - γDub[i,t] == 0)      
   else
      @constraint(m, L2dual_gDp[i in setI_D, t in setT], - πDER[i,t] + λD[B_gD[i],t] + γDlb[i,t] - γDub[i,t] == 0)
   end

   @constraint(m, L2dual_gDq[i in setI_D, t in setT], + μD[B_gD[i],t] + δDlb[i,t] - δDub[i,t] == 0)

   @constraint(m, L2dual_vD[n in setN_D, t in setT], sum(βD[l,t] for l in busesD[n].outline) - sum(βD[l,t] for l in busesD[n].inline)
      + ϑDlb[n,t] - ϑDub[n,t] == 0)
   @constraint(m, L2dual_fDp[l in setL_D, t in setT], 
      1e6(
      - sum(λD[n,t] for n in linesD[l].fbus) + sum(λD[n,t] for n in linesD[l].tbus)
      - 2*linesD[l].r * βD[l,t]
      ) == 0)
   @constraint(m, L2dual_fDq[l in setL_D, t in setT], 
      1e6(
      - sum(μD[n,t] for n in linesD[l].fbus) + sum(μD[n,t] for n in linesD[l].tbus)
      - 2*linesD[l].x * βD[l,t]
      )
      == 0)   
   @constraint(m, Utility_SD_split[t in setT], Obj_D_dual[t] - Obj_D_prim[t] <= SD_tol)      
   @constraint(m, Oprim[t in setT], Obj_W_prim[t] == 
      - sum( (Coperg[i][2] - τe) * gTp[i,t] for i in intersect(setI_T, setI_R))
      - sum( (Coperg[i][2] + CCT * H[i] ) * gTp[i,t] for i in intersect(setI_T, setI_C))    
      )      
   @constraint(m, DCPF[l in setL_T, t in setT], fTp[l,t] == (1/linesT[l].x) * (θT[linesT[l].fbus,t] - θT[linesT[l].tbus,t]))
   @constraint(m, slack[t in setT], θT[1,t] == 0)
   @constraint(m, NodeBalance_inter[n in [busT_gep], t in setT],
      sum(gTp[i,t] for i in B_gnT[n]) + sum(fTp[l,t] for l in busesT[n].inline) - sum(fTp[l,t] for l in busesT[n].outline)
      == DTp[n] * loadcurveT[n,t] + fDp0[t] - sum(gTp[i,t] for i in setI_D))
   @constraint(m, NodeBalance[n in setdiff(setN_T, busT_gep), t in setT],
      sum(gTp[i,t] for i in B_gnT[n]) + sum(fTp[l,t] for l in busesT[n].inline) - sum(fTp[l,t] for l in busesT[n].outline)
      == DTp[n] * loadcurveT[n,t])   
   @constraint(m, gpLB[i in setI_T, t in setT], gTp[i,t] >= 0)
   @constraint(m, gpUB[i in setI_T, t in setT], gTp[i,t] <= gTpo[i,t])
   @constraint(m, LineCapLB[l in setL_T, t in setT], fTp[l,t] >= - linesT[l].u)
   @constraint(m, LineCapUB[l in setL_T, t in setT], fTp[l,t] <= + linesT[l].u)
   @constraint(m, Odual[t in setT], Obj_W_dual[t] == 
      sum(gTpo[i,t]*γTub[i,t] - 0 * γTlb[i,t] for i in setI_T)
      - sum(DTp[n] * loadcurveT[n,t] * λT[n,t] for n in setdiff(setN_T, busT_gep))
      - (DTp[busT_gep] * loadcurveT[busT_gep,t] + fDp0[t] )* λT[busT_gep,t]   
      + sum(linesT[l].u * δTub[l,t] + linesT[l].u * δTlb[l,t] for l in setL_T)
      )
   @constraint(m, const_gTp_R[i in setdiff(intersect(setI_T, setI_R), setI_DER), t in setT], - γTlb[i,t] + γTub[i,t] - λT[B_g[i],t] == - (Coperg[i][2] - τe))
   @constraint(m, const_gTp_C[i in setdiff(intersect(setI_T, setI_C), setI_DER), t in setT], - γTlb[i,t] + γTub[i,t] - λT[B_g[i],t] == - (Coperg[i][2] + CCT * H[i] ))
   @constraint(m, const_gTp_DER_R[i in intersect(setI_DER, setI_R), t in setT], - γTlb[i,t] + γTub[i,t] - λT[busT_gep,t] == - (Coperg[i][2] - τe))
   @constraint(m, const_gTp_DER_C[i in intersect(setI_DER, setI_C), t in setT], - γTlb[i,t] + γTub[i,t] - λT[busT_gep,t] == - (Coperg[i][2] + CCT * H[i] ))
   @constraint(m, const_fTp[l in setL_T, t in setT], - δTlb[l,t] + δTub[l,t] + ξT[l,t] - λT[linesT[l].tbus,t] + λT[linesT[l].fbus,t] == 0)
   @constraint(m, const_θT[n in setN_T, t in setT], - sum(ξT[l,t]/linesT[l].x for l in busesT[n].outline) + sum(ξT[l,t]/linesT[l].x for l in busesT[n].inline) == 0)
   @constraint(m, WM_SD_split[t in setT], Obj_W_dual[t] - Obj_W_prim[t] <= SD_tol)
   @constraint(m, GenCostT[i in setI_T, t in setT], OCgenT[i,t] == Coperg[i][2]*gTp[i,t])
   optimize!(m)
   return m
end
