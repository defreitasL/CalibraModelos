# Calibra MD


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

hb, tp, thetb, hs =  brk_par[:,10], brk_par[:,8], brk_par[:,12], brk_par[:,7]

hb = hs

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

Hberm = 1.1

## parametros de calibración

Yerr = 5

fechas = DateTime.(YY, MM, DD, HH)
ticks = range(t_obs[1],t_obs[end],step = Year(1))

ii =  t_obs .< fechas[end]

t_obs, Y_obs_nt, Y_obs = t_obs[ii], Y_obs[ii], Y_obs_nt[ii]



w = HM.wMOORE(D50)

Omega = hb ./ (w .* tp)

idx_obs = zeros(length(t_obs))
for i in eachindex(t_obs)
    idx_obs[i] = findall((x)-> x == t_obs[i], fechas)[1]
end
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # Calculo # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Metrica a utilizar para calibrar:
MetObj = "MSS" # [Pearson, RMSE, MSS, BSS]

flagP = 4

function Calibra_MD(Χ)
    Ymd, _ = HM.MILLER_DEAN_CSonly(hb,hb./.78,sl,exp(Χ[1]),dt,D50,Hberm, exp(Χ[2]), exp(Χ[3]), Χ[4], flagP, Omega)
    YYsl = Ymd[Int.(idx_obs)]
    if MetObj == "Pearson"
        return sum((YYsl.-mean(YYsl)).*(Y_obs_nt .- mean(Y_obs_nt)))/(std(YYsl)*std(Y_obs_nt)*length(YYsl))
    elseif MetObj == "RMSE"
        return abs(sqrt(mean((YYsl .- Y_obs_nt).^2)))
    elseif MetObj == "MSS"
        return sum((YYsl .- Y_obs_nt).^2)/length(YYsl)/(var(YYsl)+var(Y_obs_nt)+(mean(YYsl)-mean(Y_obs_nt))^2)
    elseif MetObj == "BSS"
        return mean((YYsl .- Y_obs_nt .- Yerr).^2)/(Yerr).^2
    end
end

# bl = [200; 1e-5; 200]
# bu = [420; 2e-4; 720]
# P0 = [Y0; kero; kacr; YiMD; flagP]
P0 = [log(64.6031513741674);  log(0.0018894721630657854); log(0.0004669375084864312); 61.979935787939446]
mag = [P0[1] .* 1e-1; P0[2] .* 1e-1; P0[3] .* 1e-1; 4]
ngen = 10000
npop = 20
npar = length(P0)

calibr, MetVal = HM.sce_ua2(Calibra_MD, P0, ngen, npop, npar, mag)

MetVal = 1 - MetVal

if flagP ==1 println("k = k ") end
if flagP == 2 println("k = k .* hb^2") end
if flagP == 3 println("k = k .* hb^3") end
if flagP == 4 println("k = k .* Omega") end
println("----------------------------------")
println("Y0 = "*string(exp(calibr[1]))*" m")
println("kero = "*string(exp(calibr[2])))
println("kacr = "*string(exp(calibr[3])))
println("YiMD = "*string(calibr[4])*" m")
println("flagP = "*string(flagP))
println("$MetObj =  $MetVal")

Ymd, _ = HM.MILLER_DEAN_CSonly(hb,hb./.78,sl,exp(calibr[1]),dt,D50,Hberm, exp(calibr[2]), exp(calibr[3]), calibr[4], flagP, Omega)

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
    # ylim = (15, 80),
    title="Calibration of Miler and Dean model using SCE-UA algorithm | flagP = $flagP",
    titlefont=font(12, "bold")
)
plot!(
    fechas,
    Ymd,
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

savefig("MD_flagP = $flagP .png")
