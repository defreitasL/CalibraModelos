# Calibra SF


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using DataFrames
using Statistics
using DelimitedFiles
using Plots
using LaTeXStrings
using Dates
using MAT
using JLD2
using FileIO
using SciPy


# # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # Leyendo datos # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # 
wkDIR=pwd()
data_path = "C:\\Users\\freitasl\\Documents\\MATLAB\\CoastalSediments2023\\data"

include(wkDIR*"\\Modulos\\Biblioteca.jl")

data = readdlm(data_path*"\\1_PROPAGACION_ROTURA_Jara15.4_rotura06_datevec_angulo_Hcr.dat", ' ', Any, '\n'; skipblanks=true)
idx = data .!= ""
brk_par = zeros(length(idx[:,1]), sum(idx[1,:]))

for i in eachindex(idx[:,1])
    brk_par[i,:] = data[i, idx[i,:]]
end

YY, MM, DD, HH = brk_par[:,1], brk_par[:,2], brk_par[:,3], brk_par[:,4]

hb, tp, thetb =  brk_par[:,10], brk_par[:,8], brk_par[:,12]

hb[hb .< 0.05] .= 0.05

##
# ii = 1:length(idx[:,1])
# map((x) -> brk_par[x,:] = data[x, idx[x,:]], ii)
##
data = readdlm(data_path*"\\Porsmilin_tide.txt", ' ', Any, '\n'; skipblanks=true)
idx = data .!= ""
brk_par = zeros(length(idx[:,1]), sum(idx[1,:]))

for i in eachindex(idx[:,1])
    brk_par[i,:] = data[i, idx[i,:]]
end

# YY_t, MM_t, DD_t, HH_t, mm_t, tide = brk_par[:,1], brk_par[:,2], brk_par[:,3], brk_par[:,4], brk_par[:,5], brk_par[:,7]

tide = brk_par[:,7]
idx = 1:6:length(tide); tide = tide[idx]; sl = tide[1:end-1]

## Leyendo mediciones

Y0 = matread(data_path*"\\data3_4m.mat")
Y0 = Y0["data3_4m"]
t_obs = rata2datetime.(Y0["time"])
Y_obs = Y0["S3_4"]
Y_obs_nt = Y0["NEW_S3_4"]



JA15 = matread(data_path*"\\JA15_results.mat")
JA15 = JA15["JA15"]
t_JA = range(rata2datetime(JA15["timer_"][1]),rata2datetime(JA15["timer_"][end-23]),step = Hour(1))
Y_JA = JA15["xr_"][1:length(t_JA)]


# # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # Preproceso # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # 

## parametros de entrada 
dt = 1
D50 = .45e-3

Hberm = 1

## parametros de calibración

Yerr = 5

fechas = DateTime.(YY, MM, DD, HH)
ticks = range(t_obs[1],t_obs[end],step = Year(1))

ii =  t_obs .< fechas[end]

t_obs, Y_obs_nt, Y_obs = t_obs[ii], Y_obs[ii], Y_obs_nt[ii]



w = HM.wMOORE(D50)

YiSF=-10
YiMD=-20

Omega = hb ./ (w .* tp)

idx_obs = zeros(length(t_obs))
for i in eachindex(t_obs)
    idx_obs[i] = findall((x)-> x == t_obs[i], fechas)[1]
end
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # Calculo # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #



function Calibra_SF(P)
    Y, _, _, _, _, _ = HM.ShoreFor(0,tp,hb,hb./.78,D50,Omega,dt,P[1], P[2], P[3], P[4], P[5], P[6])
    YYsl = Y[Int.(idx_obs)]
    return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
end

mag = [15, 1e-3, 15, 15, 4, 0.2]
P0 = [174.57; 0.00016888956571802364; 214.558; 323.255863; -1.67025550; 0.5]
ngen = 1000
npop = 50
npar = length(P0)

calibr = CAL.sce_ua2(Calibra_SF, P0, ngen, npop, npar, mag)

println("phi = "*string(calibr[1])*" days")
println("c = "*string(calibr[2]))
println("D = "*string(calibr[3])*" days")
println("Dr = "*string(calibr[4])*" days")
println("YiSF = "*string(calibr[5])*" m")
println("k = "*string(calibr[6]))
