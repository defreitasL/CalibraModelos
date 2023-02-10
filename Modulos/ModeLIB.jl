module HM

using Polynomials
using LinearAlgebra
using SciPy
using Statistics
using ImageFiltering
using NumericalIntegration
using Shuffle


function ADEAN(D50)
    #    ###########################################################################    
    #    # Dean parameter; D50 in meters   
    #    ###########################################################################    
    A=0.51.*wMOORE(D50).^0.44
    
    return A
end

function BruunRule(hc,D50,Hberm,slr)
    #    ###########################################################################    
    #    # Bruun Rule
    #    # INPUT:
    #    # hc:     depth of closure
    #    # D50:      Mean sediment grain size (m)
    #    # Hberm:    Berm Height [m]
    #    # slr:      Expected Sea Level Rise [m]
    #    # OUTPUT:
    #    # r:        expected progradation/recession [m]
    #
    #    ###########################################################################    
    Wc = wast(hc,D50)
    
    r=slr.*Wc./(Hberm.+hc)
    
    return r
end

function depthOfClosure(Hs12,Ts12)
    #    ###########################################################################    
    #    # Closure depth, Birkemeier[1985]
    #    # Hs12:     Significant wave height exceed 12 hours in a year.
    #    # Ts12:     Significant wave period exceed 12 hours in a year.
    #    # hs_ecdf = ECDF[hsDOWserie]
    #    # f_hs = interpolate.interp1d[hs_ecdf.y,hs_ecdf.x]
    #    # Hs12 = f_hs[1-12./365./24.]
    #    # Ts12 = 5.7*np.sqrt(Hs12)
    #    ###########################################################################
            
    dc = 1.75.*Hs12 - 57.9.*(Hs12.^2. /(9.81.*Ts12.^2))
        
    return dc
end

function Hs12Calc(Hs,Tp)

    #     ###########################################################################    
    #     # Significant Wave Height exceed 12 hours a year
    #     #
    #     # INPUT:
    #     # Hs:     Significant wave height.
    #     # Tp:     Wave Peak period.
    #     #
    #     # OUTPUT:
    #     # Hs12:     Significant wave height exceed 12 hours in a year.
    #     # Ts12:     Significant wave period exceed 12 hours in a year.
    #     ###########################################################################   
    
    Hs12=zeros(size(Hs))
    Ts12=zeros(size(Tp))
    for i=1:size(Hs,2)
        Hs12calc=prctile[Hs[:,i],((365*24-12)/(365*24))*100]
        buscHS12=Hs[:,i].>=(Hs12calc-0.1) & Hs[:,i].<=(Hs12calc+0.1)
        f, xi = ksdensity[Tp[buscHS12,i]]
        _ , ii = maximum(f)
        Ts12[:,i]=xi[ii]
        Hs12[:,i]=Hs12calc
    end
    return Hs12,Ts12
end

function hunt(T,d)
    #CALCULO DE LA LONGITUD DE ONDA POR EL METODO DE HUNT
    
    ## constantes
    
    g=9.81; #[m/s^2]
    
    ## Calculos
    
    L0=g.*T.^2. /(2*pi)
    
    G=2. *pi.*(d./L0)
    
    # p=Polynomial([.067,.0864,.4622,.6522,1])
    p=Polynomial([1.0,.6522,.4622,.0864,.067])
    
    F = G .+ 1.0 ./p.(G)
    
    L=T.*(g.*d./F).^.5
    
    return L
end

function RelDisp(h, T)
    g=9.81
    
    # [h,T]=meshgrid(h,T)
    
    
    L=hunt(T,h)
    
    # Li=hunt[T,h]
    # error =1
    # while error.>1E-6
    #     L=g.*T.^2./(2*pi).*tanh(2*pi.*h./Li)
    #     error = sum(sum(abs(L-Li)))
    #     Li=L
    # end
    
    C=g.*T./(2*pi).*tanh.(2*pi.*h./L)
    
    
    
    
    
    
    return L,C
end

function RU2_Stockdon2006(slope,hs0,tp)
    #    ###########################################################################    
    #    # Run up 2# STOCKDON 2006
    #    #
    #    # INPUT:
    #    # slope:  Beach Slope in swash zone H:V. if slope is V:H all the slope terms multiply
    #    # hs0:     Significant wave height in deep water.
    #    # tp:     Peak period.
    #    #
    #    # OUTPUT:
    #    # runup2:   Run-up exceed 2#.
    
    g=9.81
    hs0 = 0.5.*(hs0[2:end]+hs0[1:end-1])
    tp = 0.5.*(tp[2:end]+tp[1:end-1])
    L0 = g.*tp.*tp./2.0/pi
    slope = 1.0./slope; # # V:H
    setup = 0.35.*slope.*(hs0.*L0).^(1.0/2.)
    infgr = (hs0.*L0.*(0.563.*slope.*slope.+0.004)).^(1.0./2.)./2.
    runup2 = 1.1 .* (setup.+infgr) # # eq 19 Stockdon 2006
    return runup2
end

function wast(hb,D50)

    #    ###########################################################################    
    #    # Width of the active surf zone
    #    # hb:   depth of closure
    #    # D50:  mean sediment grain size (m)
    #    ###########################################################################    
    #   #    hb = hb+CM ????? see why is introducing the tidal range in the width of the active surf zone    
    
    wsf=(hb./ADEAN(D50)).^(3.0/2.)
    
    return wsf
end

function wMOORE(D50)
    #    ###########################################################################    
    #    # Fall velocity; D50 in meters. Moore 1982
    #    ###########################################################################    
    
    ws = zeros(size(D50))
    for i in eachindex(D50)
        if D50[i].<=0.1*1e-3
            ws[i] = 1.1*1e6*D50[i].^2
        elseif D50[i].>0.1*1e-3 && D50[i].<1.0*1e-3
            ws[i] = 273.0.*D50[i].^1.1
        elseif D50[i].>1*1e-3
            ws[i] = 4.36.*D50[i].^0.5
        end
    end
    
    return ws
end

function Y09_MOD(hs,HSB,DY,DT)

    tau = DT .* (hs ./ HSB).^-1

    Yeq = -DY .* (hs.^2 - HSB.^2)./HSB.^2

    return Yeq, tau

end

function ShoreFor(OmegaEQ,tp,hb,depthb,D50,Omega,dt,phi = 0, c = 0, D = 0, Dr = 0, Sini = 0, k = 0.5)
    rho = 1025.
    g = 9.81
    if length(size(hb)) != 1
        
        hbr = .5 .*(hb[2:end]+hb[1:end-1])
        # hbr[1] = hbr[2]
        # hbr[end] = hbr[end-1]
        depthbr = .5 .*(depthb[2:end]+depthb[1:end-1])
        # depthbr[1] = depthbr[2]
        # depthbr[end] = depthbr[end-1]
        tpr = .5 .*(tp[2:end]+tp[1:end-1])
        # tpr[1] = tpr[2]
        # tpr[end] = tpr[end-1]
    end
    if size(Omega,2) == 3
        P = 1 ./16 .*ro .* g .* hbr.^2 .* (g.*depthbr).^.5
        ws = wMOORE(D50)
        OmegaN = hbr ./ (ws .* tpr)
        F = P.^.5 .* (OmegaEQ - OmegaN)./dpO

        return F, OmegaN, OmegaEQ

    elseif length(size(hb)) != 1
        P = 1 ./16 .*ro .* g .* hbr.^2 .* (g.*depthbr).^.5
        ws = wMOORE(D50)
        Omega[:,end] = hbr ./ (ws .* tpr)
        phi = size(Omega,2)/2/24*dt
        phivec = 1:dt/24:2*phi+(24-dt)/24
        OmegaEQ = zeros(size(OmegaEQ))
        for i in eachindex(Omega[:,1])
            OmegaEQ[i] = sum(Omega[i,:].*phivec)./sum(phivec)
        end
        F = P.^.5 .* (OmegaEQ - Omega[:,end])./dpO

        return F, Omega[:,end], OmegaEQ
    else
        # ddt = dt/24
        # phi = phi*24
        P = 1 ./16 .*rho .* g .* hb.^2 .* (g.*depthb).^.5
        ii = 1:dt:D*24 #-dt*24
        phivecP = 10 .^(-abs.(ii)./(phi*24))
        # phivecP = reverse(10 .^(-abs.(ii)./phi))
        IDX = length(phivecP)
        # println(phivecP)

        
        
        
        # OmegaEQ = zeros(length(Omega)-IDX)
        
        # F = zeros(length(Omega)-IDX)
        # dS_dt = zeros(length(Omega)-IDX)
        # S = zeros(length(Omega)-IDX)
        # PHISUM = sum(phivecP)
        # for i in eachindex(OmegaEQ)
        #     @fastmath OmegaEQ[i] = sum(Omega[i:i-1+IDX].*phivecP)./PHISUM
        #     # @fastmath F[i] = P[IDX+i-1].^.5 .* (OmegaEQ[i] - Omega[IDX+i-1])./dpO
        #     # if F[i] > 0
        #     #     F_ero = 0
        #     #     F_acr = F[i]
        #     # else
        #     #     F_ero = F[i]
        #     #     F_acr = 0
        #     # end
        #     # @fastmath dS_dt[i] = c_ero .* F_ero  + c_acr .* F_acr
        #     # if i ==1
        #     #     @fastmath S[i] = dt * dS_dt[i] + Sini
        #     # else
        #     #     @fastmath S[i] = dt * dS_dt[i] + S[i-1]
        #     # end
        # end
        # F = P[IDX+1:end].^.5 .* (OmegaEQ .- Omega[IDX+1:end])./std(OmegaEQ)
        

        # println("D = " *string(D))
        # println("phi = " *string(phi))
        # println("c = " *string(c))
        # reflect(centered(vent))
        OmegaAUX = Omega .- mean(Omega)
        phivecP = [zeros(IDX-1); phivecP]
        vent = reflect(centered(phivecP./sum(phivecP)))
        OmegaEQ = imfilter(OmegaAUX,vent, Fill(0))
        # println(size(DSP.conv(OmegaAUX,window)))
        # OmegaEQ = OmegaEQ[IDX:end] .+ mean(Omega)
        # F = P[IDX:end].^.5 .* (OmegaEQ .- Omega[IDX:end])./dpO
        OmegaEQ = OmegaEQ .+ mean(Omega)
        F = P.^k .* (OmegaEQ .- Omega)./std(OmegaEQ)
        S = zeros(length(Omega))
        # S = zeros(length(Omega)-IDX)
        
        # Fdt = SciPy.signal.detrend(F) .+ mean(F)

        N = length(1:dt:Dr*24)
        rero = F .< 0
        racr = F .> 0
        

        # r= abs(sum((Fdt).* racr)/sum((Fdt).* rero))
        
        # r= abs(sum(F[racr].- SciPy.signal.detrend(F[racr]).+ mean(F[racr]))/sum(F[rero].- SciPy.signal.detrend(F[rero]).+ mean(F[rero])))
        # println("r = "*string(r))

        # r = abs(sum(F[racr])/sum(F[rero]))
        # @fastmath dS_dt = c.* r .* F .*rero  + c .* F .* racr
        
        S[1] = Sini
        
        # for i in eachindex(S[2:end])
        #     S[i+1] = dt * dS_dt[i+1] + S[i]
        # end

        # r = abs.(Facr./Fero)
        r = zeros(length(Omega))
        r[1:N] .= abs(sum(F[1:N].*racr[1:N])/sum(F[1:N].*rero[1:N]))
        Fcacr = cumulative_sum(F.*racr, N)
        Fcero = cumulative_sum(F.*rero, N)
        
        r[N:end] .= abs.(Fcacr./Fcero)
        # mmF = moving_average(F, N)
        # Fcacr = cumsum(F .* racr)
        # Fcero = cumsum(F .* rero)
        # for i in eachindex(S[2:end])
        #     # if i > N
        #     #     r[i+1] = abs(sum((F[i-N:i+1]) .* racr[i-N:i+1])/sum((F[i-N:i+1])  .*rero[i-N:i+1])) #i-N
        #     #     # r[i+1] = abs(Fcacr[i+1]/Fcero[i+1])
        #     # end
        #     S[i+1] = dt *(c * 0.5 * (r[i+1] .* rero[i+1] .* F[i+1] + racr[i+1] .* F[i+1] + r[i] .* rero[i] .* F[i] + racr[i] .* F[i]) )+ S[i]
        #     # S[i+1] = dt *(c * 0.5 * (r .* rero[i+1] .* F[i+1] + racr[i+1] .* F[i+1] + r .* rero[i] .* F[i] + racr[i] .* F[i]) )+ S[i]
        # end
        

        dt_c_half = 0.5 * dt * c
        r_rero_F = r[2:end] .* rero[2:end] .* F[2:end]
        racr_F = racr[2:end] .* F[2:end]
        r_rero_F_prev = r[1:end-1] .* rero[1:end-1] .* F[1:end-1]
        racr_F_prev = racr[1:end-1] .* F[1:end-1]
        S[2:end] = dt_c_half .* cumsum(r_rero_F .+ racr_F .+ r_rero_F_prev .+ racr_F_prev) .+ S[1]

        dS_dt = zeros(length(Omega)-IDX)

        return S, F, P, OmegaEQ, dS_dt, length(Omega) - length(S)+1
    end
    
end

function MILLER_DEAN_CSonly(hb,depthb,sl,Y0,dt,D50,Hberm, kero, kacr, Yi, flagP = 1, Omega = 0)

    if flagP == 1
        kero = repeat(kero, length(hb))
        kacr = repeat(kacr, length(hb))
    elseif flagP == 2
        kero = kero .* hb .^2
        kacr = kacr .* hb .^2
    elseif flagP == 3
        kero = kero .* hb .^3
        kacr = kacr .* hb .^3
    elseif flagP == 4
        kero = kero .* Omega
        kacr = kacr .* Omega
    end
    yeq = zeros(size(hb))
    Y = zeros(size(hb))

    for i in eachindex(sl)
        wl = 0.106.*hb[i].+sl[i] # # total water level
        Wast = wast(depthb[i],D50)
        yeq[i] = Y0 .- Wast.*(wl)./(Hberm.+depthb[i])
        if i ==1
            r = yeq[i].-Y[i] > 0
            k = kacr[i] * r + kero[i] * !r
            # k = k .* hb[i]^2
            Y[i] = Yi
        else
            r = yeq[i].-Y[i-1] > 0
            k = kacr[i] * r + kero[i] * !r
            # k = k .* hb[i]^2
            Y[i] = dt * k * (yeq[i] - Y[i-1]) + Y[i-1]
        end
    end        
    return Y, yeq


end

#################################################################

function ShoreForCAL(dpO,hb,depthb,Omega,dt,phi = 0, c = 0, r = 1, Sini = 0)
    ro = 1025.
    g = 9.81
    # ddt = dt/24
    phi = phi*24
    HD = 2*phi
    P = 1 ./16 .*ro .* g .* hb.^2 .* (g.*depthb).^.5
    ii = 0:HD-1
    phivecP = 10 .^(-abs.(ii)./phi)
    # phivecP = reverse(10 .^(-abs.(ii)./phi))
    IDX = length(phivecP)

    OmegaAUX = Omega .- mean(Omega)
    phivecP = [zeros(IDX-1); phivecP]
    vent = phivecP./sum(phivecP)
    S = zeros(length(Omega))
    # OmegaEQ = DSP._conv_similar(OmegaAUX,window, size(OmegaAUX))
    OmegaEQ = imfilter(OmegaAUX,reflect(centered(vent)), Fill(0))
    # println(size(DSP.conv(OmegaAUX,window)))
    OmegaEQ = OmegaEQ .+ mean(Omega)
    F = P.^.5 .* (OmegaEQ .- Omega)./dpO
    # F = P.^.5 .* (OmegaEQ .- Omega)./dpO

    rero = F .< 0
    racr = F .> 0

    @fastmath dS_dt = c.* r .* F .*rero  + c .* F .* racr
    S[1] = dt * dS_dt[1] + Sini
    for i in eachindex(S[2:end])
        S[i+1] = dt * dS_dt[i+1] + S[i]
    end
    
    return S

    
end

function MILLER_DEAN_CSonlyCAL(hb,depthb,sl,Y0,dt,D50,Hberm, kacr, kero)

        
    yeq = zeros(size(hb))
    Y = zeros(size(hb))

    for i in eachindex(sl)
        wl = 0.106.*hb[i].+sl[i]; # # total water level
        Wast = wast(depthb[i],D50)
        yeq[i] = Y0 .- Wast.*(wl)./(Hberm.+depthb[i])
        
        if i ==1
            r = yeq[i].-Y[i] > 0
            k = kacr * r + kero * !r
            Y[i] = dt * k * (yeq[i] - Y0)
        else
            r = yeq[i].-Y[i-1] > 0
            k = kacr * r + kero * !r
            Y[i] = dt * k * (yeq[i] - Y[i-1]) + Y[i-1]
        end
    end        
    return Y, yeq


end

########################################################

function moving_average(data, window_size)

    result = zeros(length(data))
    result[1:length(data) - window_size] .= mean(data[1:length(data) - window_size])
    for i = 1:(length(data) - window_size + 1)
        result[i] = sum(data[i:(i + window_size - 1)]) / window_size
    end

    return result
end

function cumulative_sum(input_array, window_size)
    n = length(input_array)
    cum_sum = cumsum(input_array)
    # output = zeros(n - window_size + 1)
    # for i in 1:(n - window_size + 1)
    #     output[i] = cum_sum[i + window_size - 1] - (i - 1 >= 1 ? cum_sum[i - 1] : 0)
    # end
    return cum_sum[window_size:end] .- [0; cum_sum[1:(n-window_size)]]
end





end # module
# Sure! Here's an explanation of each input variable in the SCE-UA algorithm in Julia:

# f: This is the function that you want to optimize. It should take a single input argument x (of length npar) and return a scalar value.

# x0: This is the initial guess for the solution to the optimization problem. It should be a 1-dimensional array of length npar.

# ngen: This is the number of generations (i.e., iterations) of the SCE-UA algorithm.

# npop: This is the number of individuals (candidate solutions) in the population. It should be an even number, as the SCE-UA algorithm requires that the population size be divisible by 2.

# npar: This is the number of parameters in the optimization problem. It should be the length of x0 and the number of elements in each row of the pop array.