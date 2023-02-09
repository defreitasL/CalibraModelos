# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# set JULIA_NUM_THREADS = 6

using DataFrames
using Statistics
# using KrylovKit ### linsolve
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

level = 0
Yerr = 5

figname = string(level)*"_Calibration"


#ShoreFor
phi = 157.28
c = 1.607e-4
D = 228.96
Dr = 317.71
YiSF=-0.947

# r = SF["r"][findmax(SF["MS"])[2][3]]
# M&D
Y0 = 25
kero = 0.00558106796640481
kacr = 0.0003214400492900498

YiMD = -15.23444928731412
#


fechas = DateTime.(YY, MM, DD, HH)
ticks = range(t_obs[1],t_obs[end],step = Year(1))

ii =  t_obs .< fechas[end]

t_obs, Y_obs_nt, Y_obs = t_obs[ii], Y_obs[ii], Y_obs_nt[ii]



w = HM.wMOORE(D50)



# Yi = Y_obs[1]
# id = findall((x)-> x == t_obs[1], fechas)[1]

# tp, hb, sl, fechas = tp[id:end], hb[id:end], sl[id:end], fechas[id:end]


Omega = hb ./ (w .* tp)

idx_obs = zeros(length(t_obs))
for i in eachindex(t_obs)
    idx_obs[i] = findall((x)-> x == t_obs[i], fechas)[1]
end
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # Calculo # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #




Y, _, _, _, _, iidd = HM.ShoreFor(0,tp,hb,hb./.78,D50,Omega,dt,phi, c, D, Dr, YiSF)
YYsl = Y[Int.(idx_obs)]
RP1 = sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))
MS1 = 1 - sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
BSS1 = 1 - mean((YYsl .- Y_obs .- Yerr).^2)/(Yerr).^2
RMSE1 = sqrt(mean((YYsl .- Y_obs).^2))
println("ShoreFor - OK")
println("RMSE = " , RMSE1)
println("RP = " , RP1)
println("MS = " , MS1)
println("BSS = " , BSS1)



Ymd, _ = HM.MILLER_DEAN_CSonly(hb,hb./.78,sl,Y0,dt,D50,Hberm, kero, kacr, YiMD)
YYsl = Ymd[Int.(idx_obs)]
RP2 = sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))
RMSE2 = sqrt(mean((YYsl .- Y_obs).^2))
MS2 = 1 - sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
BSS2 = 1 - mean((YYsl .- Y_obs .- Yerr).^2)/(Yerr).^2
println("Miller y Dean - OK")
println("RMSE = " , RMSE2)
println("RP = " , RP2)
println("MS = " , MS2)
println("BSS = " , BSS2)

idx_obs = zeros(length(t_obs))
for i in eachindex(t_obs)
    idx_obs[i] = findall((x)-> x == t_obs[i], t_JA)[1]
end
YYsl = Y_JA[Int.(idx_obs)]
RP3 = sum((YYsl.-mean(YYsl)).*(Y_obs_nt .- mean(Y_obs_nt)))/(std(YYsl)*std(Y_obs_nt)*length(YYsl))
RMSE3 = sqrt(mean((YYsl .- Y_obs_nt).^2))
MS3 = 1 - sum((YYsl .- Y_obs_nt).^2)/length(YYsl)/(var(YYsl)+var(Y_obs_nt)+(mean(YYsl)-mean(Y_obs_nt))^2)
BSS3 = 1 - mean((YYsl .- Y_obs_nt .- Yerr).^2)/(Yerr).^2
println("Jara")
println("RMSE = " , RMSE3)
println("RP = " , RP3)
println("MS = " , MS3)
println("BSS = " , BSS3)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # Plots # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #


scatter(
    t_obs,
    Y_obs_nt,
    label="Datos", 
    tickfont = (10,"Computer Modern"),
    fontfamily = "Computer Modern",
    ylabel = "Shoreline Position (m)",
    labelfontsize = 12,
    left_margin = 20Plots.mm,
    right_margin = 8Plots.mm,
    bottom_margin = 5Plots.mm,
    size = [2000, 400],
    markershape = :circle,
    markeralpha = 0.8,
    markercolor = :gray,
    markerstrokewidth = 1,
    markerstrokealpha = 1,
    markerstrokecolor = :black,
    legend=:topleft,
    xlim = (DateTime(2004),DateTime(2020)),
    ylim = (15, 80)
)
plot!(
    t_JA,
    Y_JA,
    label="Jara et al.",
    lc=:black,
    lw=1
)
plot!(
    fechas[iidd:end],
    Y.+ mean(Y_obs_nt),
    label="ShoreFor",
    lc=:red,
    lw=1,
)
plot!(
    fechas,
    Ymd .+ mean(Y_obs_nt),
    label="Miller & Dean",
    lc=:blue,
    lw=1
)

plot!(twinx(),
    fechas,
    hb.^2 ./maximum(hb.^2.),
    lc = :magenta,
    lw = 0.8,
    label = "Wave Energy",
    legend=:topright,
    ylim = (0, 4),
    xticks = ([0:0.25:1], [0:0.25:1]))
plot!(xticks = (ticks, Dates.format.(ticks, "yyyy")))
# annotate!(fechas[24*180], 0, ("RMSE = ", 7, :black), :black, )

savefig(figname * ".png")

