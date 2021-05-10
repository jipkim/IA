using Pkg, DataFrames, Distributions, DelimitedFiles, CSV
using JuMP, Gurobi, KNITRO
Solver = "KNITRO"
setT = 1:24

κ = 0.7
testsystemT = "NYISO_v2_70_rCoal"
testsystemD = "Manhattan"

σ0_pv = 20e-2
σ0_wind = 20e-2/2
υ0 = 0.0
UserGap = 1e-6
D_ratio = 0.2
D_scale = 1
D_new = 1
busT_gep = 10
busD_gep = 4
SD_tol = 0
DCRF = round(0.0002198427046320309 * length(setT) / 24, digits=5)
ScaleObj = 1

semi_cut_MP = 0
semi_cut_SP = 0
round_MP = 1
round_SP = 0

policy_switch = 1
IA_switch = 1
IA_LMP = 0 
WMP = 0 
CTI = 0 
ED =  50
γCT = 1
WHA = 0


η = 0.03

include("StructType_v2.jl")
bound = bounds(0, 200, 0, 200, 0, 200, 0, 200, 0, 0)

include("MPini.jl")
include("MP.jl")   
include("SP.jl")

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

spi_ratio = 0.05
π0 = 50
s0 = spi_ratio * π0
d0 = DDp * loadcurveD'
N = 0.25
M = N * d0 .+ π0 .+ s0
M_new = 11(π0+s0)
N_new = 10(π0+s0)ones(size(d0))./d0
ST_res = s0 * ones(length(setN_T), maximum(setT)) # surcharge for transmission nodes

DTp_inter = DTp[busT_gep] * loadcurveT[busT_gep,:]
M_new = 11(π0+s0)
N_new_DTp = 10(π0+s0)ones(size(DTp_inter))./DTp_inter



time_M = @elapsed tempM = MPini(testsystemT, testsystemD, busT_gep, busD_gep, κ, setT, η, σ0_pv, σ0_wind, υ0, Solver, bound, D_ratio, ScaleObj, SD_tol, M, N, ED, γCT, ST_res, WMP, CTI)
if round_SP == 1
   time_S = @elapsed tempS = SP(testsystemT, testsystemD, busT_gep, busD_gep, κ, setT, η, σ0_pv, σ0_wind, υ0, 
      round.(value.(tempM[:πDER]), digits=5), round.(value.(tempM[:πC]), digits=5), round.(value.(tempM[:s]), digits=5),
      round.(value.(tempM[:τe]), digits=5), round.(value.(tempM[:τc]), digits=5), 
      Solver, D_ratio, ScaleObj, SD_tol, M, N, ED, γCT, ST_res, WMP, CTI)
else
   time_S = @elapsed tempS = SP(testsystemT, testsystemD, busT_gep, busD_gep, κ, setT, η, σ0_pv, σ0_wind, υ0, value.(tempM[:πDER]), value.(tempM[:πC]), value.(tempM[:s]),
      value.(tempM[:τe]), value.(tempM[:τc]), Solver, D_ratio, ScaleObj, SD_tol, M, N, ED, γCT, ST_res, WMP, CTI)
end
Niter = 100
conv = 777*ones(Niter,4)
time_iter = zeros(Niter,2)
MPbound_iter = zeros(Niter)
SPbound_iter = zeros(Niter)

πDER_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempM[:πDER])}}(undef,Niter)
πC_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempM[:πC])}}(undef,Niter)
s_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempM[:s])}}(undef,Niter)
τe_iter = 777*ones(Niter)
τc_iter = 777*ones(Niter)
ObjReg_M_iter = 777*ones(Niter)
ObjDER_M_iter = 777*ones(Niter)
ObjDM_M_iter = 777*ones(Niter)
ObjWM_M_iter = 777*ones(Niter)
ObjReg_S_iter = 777*ones(Niter)
ObjDER_S_iter = 777*ones(Niter)
ObjDM_S_iter = 777*ones(Niter)
ObjWM_S_iter = 777*ones(Niter)

gpmax_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gpmax])}}(undef,Niter)
gTpo_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gTpo])}}(undef,Niter)
gDpo_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gDpo])}}(undef,Niter)
gDqo_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gDqo])}}(undef,Niter)
gTp_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gTp])}}(undef,Niter)
gDp_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gDp])}}(undef,Niter)
gDq_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gDq])}}(undef,Niter)
λT_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:λT])}}(undef,Niter)
λD_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:λD])}}(undef,Niter)
status_iter = Array{Union{Nothing, String}}(nothing, Niter, 2)

ii = 1
time_iter[ii,1] = time_M
time_iter[ii,2] = time_S
status_iter[ii,1] = string(termination_status(tempM))
status_iter[ii,2] = string(termination_status(tempS))

πDER_iter[ii] = value.(tempM[:πDER])
πC_iter[ii] = value.(tempM[:πC])
s_iter[ii] = value.(tempM[:s])
τe_iter[ii] = value.(tempM[:τe])
τc_iter[ii] = value.(tempM[:τc])
ObjReg_M_iter[ii] = value.(tempM[:Obj_REG])
ObjDER_M_iter[ii] = value.(tempM[:Obj_DER])
ObjDM_M_iter[ii] = sum(value.(tempM[:Obj_D_prim_act]))
ObjWM_M_iter[ii] = sum(value.(tempM[:Obj_W_prim]))
ObjReg_S_iter[ii] = value.(tempS[:Obj_REG])
ObjDER_S_iter[ii] = value.(tempS[:Obj_DER])
ObjDM_S_iter[ii] = sum(value.(tempS[:Obj_D_prim_act]))
ObjWM_S_iter[ii] = sum(value.(tempS[:Obj_W_prim]))

if round_MP == 1
   gpmax_iter[ii] = abs.(round.(value.(tempS[:gpmax]), digits=5))
   gTpo_iter[ii] = abs.(round.(value.(tempS[:gTpo]), digits=5))
   gDpo_iter[ii] = abs.(round.(value.(tempS[:gDpo]), digits=5))
   gDqo_iter[ii] = abs.(round.(value.(tempS[:gDqo]), digits=5))
else
   gpmax_iter[ii] = abs(value.(tempS[:gpmax]))
   gTpo_iter[ii] = abs(value.(tempS[:gTpo]))
   gDpo_iter[ii] = abs(value.(tempS[:gDpo]))
   gDqo_iter[ii] = abs(value.(tempS[:gDqo]))
end


gTp_iter[ii] = value.(tempS[:gTp])
gDp_iter[ii] = value.(tempS[:gDp])
gDq_iter[ii] = value.(tempS[:gDq])
λT_iter[ii] = value.(tempS[:λT])
λD_iter[ii] = value.(tempS[:λD])

conv_temp = (value.(tempM[:Obj_REG]) - value.(tempS[:Obj_REG]))
global conv[ii,:] = [round.(value.(tempM[:Obj_REG]), digits = 3), round.(value.(tempS[:Obj_REG]), digits = 3),
   round.(value.(tempM[:Obj_REG]) - value.(tempS[:Obj_REG]), digits = 3), (value.(tempM[:Obj_REG]) - value.(tempS[:Obj_REG])) / value.(tempM[:Obj_REG])]
global SPbound_iter[ii] = value.(tempS[:SPbound])
global MPbound_iter[ii] = value.(tempM[:MPbound])
df = DataFrame(τc = round.(τc_iter[1:ii], digits = 3), τe = round.(τe_iter[1:ii], digits = 3),
   ObjReg_M = round.(ObjReg_M_iter[1:ii], digits=3), 
   ObjDER_M = round.(ObjDER_M_iter[1:ii], digits=3),
   ObjDM_M = round.(ObjDM_M_iter[1:ii], digits=3), 
   ObjWM_M = round.(ObjWM_M_iter[1:ii], digits=3),
   ObjReg_S = round.(ObjReg_S_iter[1:ii], digits=3), 
   ObjDER_S = round.(ObjDER_S_iter[1:ii], digits=3),
   ObjDM_S = round.(ObjDM_S_iter[1:ii], digits=3), 
   ObjWM_S = round.(ObjWM_S_iter[1:ii], digits=3)
   )

for ii in 2:Niter   
   time_M = @elapsed global tempM = MP(testsystemT, testsystemD, busT_gep, busD_gep, κ, setT, η, σ0_pv, σ0_wind, υ0, Solver, bound, D_ratio, ScaleObj, SD_tol, M, N, ED, γCT, ST_res, WMP, CTI, gpmax_iter, gTpo_iter, gDpo_iter, gDqo_iter, ii-1)
   global time_iter[ii,1] = time_M
   global status_iter[ii,1] = string(termination_status(tempM))

   global πDER_iter[ii] = value.(tempM[:πDER])
   global πC_iter[ii] = value.(tempM[:πC])
   global s_iter[ii] = value.(tempM[:s])
   global τe_iter[ii] = value.(tempM[:τe])
   global τc_iter[ii] = value.(tempM[:τc])
   global ObjReg_M_iter[ii] = value.(tempM[:Obj_REG])   
   global ObjDER_M_iter[ii] = value.(tempM[:Obj_DER])
   global ObjDM_M_iter[ii] = sum(value.(tempM[:Obj_D_prim_act]))
   global ObjWM_M_iter[ii] = sum(value.(tempM[:Obj_W_prim]))

   
   time_S = @elapsed global tempS = SP(testsystemT, testsystemD, busT_gep, busD_gep, κ, setT, η, σ0_pv, σ0_wind, υ0, value.(tempM[:πDER]), value.(tempM[:πC]), value.(tempM[:s]),value.(tempM[:τe]), value.(tempM[:τc]), Solver, D_ratio, ScaleObj, SD_tol, M, N, ED, γCT, ST_res, WMP, CTI)
      
   global time_iter[ii,2] = time_S
   global status_iter[ii,2] = string(termination_status(tempS))

   global ObjReg_S_iter[ii] = value.(tempS[:Obj_REG])
   global ObjDER_S_iter[ii] = value.(tempS[:Obj_DER])
   global ObjDM_S_iter[ii] = sum(value.(tempS[:Obj_D_prim_act]))
   global ObjWM_S_iter[ii] = sum(value.(tempS[:Obj_W_prim]))
   global gpmax_iter[ii] = abs(value.(tempS[:gpmax]))
   global gTpo_iter[ii] = abs(value.(tempS[:gTpo]))
   global gDpo_iter[ii] = abs(value.(tempS[:gDpo]))
   global gDqo_iter[ii] = abs(value.(tempS[:gDqo]))
   
   global gTp_iter[ii] = value.(tempS[:gTp])
   global gDp_iter[ii] = value.(tempS[:gDp])
   global gDq_iter[ii] = value.(tempS[:gDq])
   global λT_iter[ii] = value.(tempS[:λT])
   global λD_iter[ii] = value.(tempS[:λD])

   conv_temp = (value.(tempM[:Obj_REG]) - value.(tempS[:Obj_REG]))
   global conv[ii,:] = [round.(value.(tempM[:Obj_REG]), digits = 3), round.(value.(tempS[:Obj_REG]), digits = 3),
      round.(value.(tempM[:Obj_REG]) - value.(tempS[:Obj_REG]), digits = 3), (value.(tempM[:Obj_REG]) - value.(tempS[:Obj_REG])) / value.(tempM[:Obj_REG])]
   global SPbound_iter[ii] = value.(tempS[:SPbound])
   global MPbound_iter[ii] = value.(tempM[:MPbound])
   df = DataFrame(τc = round.(τc_iter[1:ii], digits = 3), τe = round.(τe_iter[1:ii], digits = 3),
      ObjReg_M = round.(ObjReg_M_iter[1:ii], digits=3), 
      ObjDER_M = round.(ObjDER_M_iter[1:ii], digits=3),
      ObjDM_M = round.(ObjDM_M_iter[1:ii], digits=3), 
      ObjWM_M = round.(ObjWM_M_iter[1:ii], digits=3),
      ObjReg_S = round.(ObjReg_S_iter[1:ii], digits=3), 
      ObjDER_S = round.(ObjDER_S_iter[1:ii], digits=3),
      ObjDM_S = round.(ObjDM_S_iter[1:ii], digits=3), 
      ObjWM_S = round.(ObjWM_S_iter[1:ii], digits=3)
   )   
end
