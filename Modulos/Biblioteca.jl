module HM

using Polynomials
using LinearAlgebra
using SciPy
using Statistics
using ImageFiltering
using NumericalIntegration
using Shuffle

function abs_angle_cartesian(relD, batiD)
    #     ###########################################################################    
    #     # Absolute angle in cartesian notation, angle between [180,-180], 
    #     # 0 is in EAST & positive counterclockwise.
    #     # From a relative angle from wave & bathymetry.
    #     # The same as rel_angle_cartesian[relD,-1*batiD]
    #     # INPUT:
    #     # relD:     relative wave angle between wave and bathymetry; 0 is the bathymetry & positive counterclockwise.
    #     # batiD:    bathymetry angle (normal to the shoreline) in Cartesian notation.
    #     #
    #     # OUTPUT:
    #     # waveD:    wave angle in Cartesian notation.
    #     ###########################################################################    
    
    waveD = relD + batiD
    replace!(waveD,)
    waveD[waveD.>180]=waveD[waveD.>180].-360
    waveD[waveD.<-180]=waveD[waveD.<-180].+360
        
    return waveD
end

function abs_pos(X0,Y0,phi,dn)
    #     #####################    
    #     # INPUT:
    #     #
    #     # X0 : x coordinate; origin of the transect
    #     # Y0 : y coordinate; origin of the transect
    #     # phi : transect orientation in radians
    #     # dn : position on the transect
    #     #
    #     # OUTPUT:
    #     # XN : x coordinate
    #     # YN : y coordinate
    #     #####################
    # n=size(dn,2)
    XN = X0 + dn.*cos.(phi)
    YN = Y0 + dn.*sin.(phi)
       
    #     return XN; YN
  
    return XN, YN
end

function ADEAN(D50)
    #    ###########################################################################    
    #    # Dean parameter; D50 in meters   
    #    ###########################################################################    
    A=0.51.*wMOORE(D50).^0.44
    
    return A
end

function ALST(Hb,Tb,Dirb,hb,ANGbati,K)


    # q=ALST(Hb,Tb,Dirb,hb,ANGbati,method,K)
    #    ###########################################################################    
    #    # Alongshore sediment transport
    #    #    
    #    # INPUT:
    #    # Hb:        wave height.
    #    # Tb:        wave period.
    #    # Dirb:      wave direction. Nautical convention.
    #    # hb:        depth of wave conditions.
    #    # ANGbati:   bathymetry angle; the normal of the shoreline in nautical degrees.
    #    # method:    alongshore sediment transport formula. Default is CERQ
    #    # calp:      calibration parameters
    #    # K1:        Komar calibration parameter
    #    # K:         SPM calibration parameter
    #    #    
    #    # OUTPUT:
    #    # q:      alongshore sediment transport relative to the bathymetry angle.
    #    #
    #    # DEPENDENCIAS:
    #    # rel_angle_cartesian; DIRrel
    #    ###########################################################################
    #    #Hb=hb[:,ti] Tb=tb[:,ti]" Dirb=dirb[:,ti] hb=depthb[:,ti]";ANGbati=alfas+90;method=sedtr;K=ka
    
    DIRrel = rel_angle_cartesian(nauticalDir2cartesianDir(Dirb),ANGbati)
    PerpRange = abs.(DIRrel) .< 90
    PerpRange = vec(PerpRange)
    q = zeros(size(Hb))
    q0 = zeros(size(Hb))
    # K1 = 0.39; # dimensionless, 0.39 SPM using Hs &  0.77 using Hrms Komar[1970]
    # K1 = calp.K1
    # q[PerpRange] = K[PerpRange].*(Hb[PerpRange].^2).*...
    #     cos(deg2rad[DIRrel[PerpRange]]).*sin(deg2rad[DIRrel[PerpRange]])
    rho = 1025 # #saltwater mass density SPM
    rhos = 2650 # # sand mass density SPM
    
    p = 0.4 # # porosity SPM
    gammab = Hb[PerpRange]./hb[PerpRange]
    gammab[isnan.(gammab)].=Inf
    cnts = rho.*sqrt.(9.81)./(16. .*sqrt.(gammab).*(rhos.-rho).*(1.0.-p));      
    q0[PerpRange] .= K[PerpRange].*cnts.*(Hb[PerpRange].^(5. /2.))
    q[PerpRange] .= q0[PerpRange].*sin.(2. .*deg2rad.(DIRrel[PerpRange]))
    #     elif method .== "spm":
    # #        K = 0.0098
    #         rho = 1025
    #         rhos = 2650
    #         lbda = 0.4
    # #        K = calp.K
    #         q[PerpRange] = rho/8.*K[PerpRange]/((rhos-rho)*lbda)*Hb[PerpRange]*np.sqrt(9.81*hb[PerpRange])*np.cos(np.radians[DIRrel[PerpRange]])*np.sin(np.radians[DIRrel[PerpRange]])
    #         
    # q[1] = q[2] + (q[3] - q[2])

    # if q[3]-q[2] < 0
    #     q[1] = q[2] + (q[3] - q[2])
    # else
    #     q[1] = q[2] - (q[3] - q[2])
    # end

    # q[end] = q[end-1] + (q[end-2] - q[end-1])
    # if q[end-2]-q[end-1] < 0
    #     q[end] = q[end-1] + (q[end-2] - q[end-1])
    # else
    #     q[end] = q[end-1] - (q[end-2] - q[end-1])
    # end


    q[1]=q[2]
    q[end]=q[end-1]
    # q[1]=0
    # q[end]=0
    

    return q, q0
end

function BreakingPropagation(H1,T1,DIR1,h1,ANGbati)
    # BreakingPropagation(H1,T1,DIR1,h1,ANGbati,breakType)
    #     ###########################################################################    
    #     # Propagation of waves using linear theory assuming rectilinear & parallel bathymetry
    #     #    
    #     # INPUT:
    #     # H1:        wave height.
    #     # T1:        wave period.
    #     # DIR1:      wave direction. Nautical convention.
    #     # h1:        depth of wave conditions.
    #     # ANGbati:   bathymetry angle; the normal of the shoreline. Cartesian notation
    #     # breakType: type of breaking condition. Spectral | monochromatic.
    #     #    
    #     # OUTPUT:
    #     # H2:        wave height during breaking. Wave period is assumed invariant due to linear theory
    #     # DIR2:      wave direction during breaking. Nautical convention.
    #     # h2:        depth of breaking
    #     ###########################################################################    
    
    # breaksT = ("spectral':0.45, 'mono":0.78)
    
    # Bcoef=breaksT[breakType.lower()]
        
    Bcoef=0.55
    DIRrel = rel_angle_cartesian(nauticalDir2cartesianDir(DIR1),ANGbati)
    
    
    h2l0 = H1./Bcoef; # # initial condition for breaking depth
        
    
    H2 = zeros(length(H1),1); 
    DIR2 = zeros(length(DIR1),1); 
    h2 = zeros(length(H1),1)
    
    
    H2[h2l0.>=h1] = H1[h2l0.>=h1]
    DIR2[h2l0.>=h1] = DIR1[h2l0.>=h1]
    h2[h2l0.>=h1] = h2l0[h2l0.>=h1]  # # check that the initial depth is deeper than the breaking value
    
    H2[H1.<=0.1] = H1[H1.<=0.1]
    DIR2[H1.<=0.1] = DIR1[H1.<=0.1]
    h2[H1.<=0.1] = h2l0[H1.<=0.1]  # # check that the initial depth is deeper than the breaking value
    
    propProf = (abs.(DIRrel).<=90) .&& (H1.>0.1) .&& (h2l0.<h1)
    propProf = vec(propProf)
    # #    print "init"    
    # #    print len[propProf]
    # #    print "end"
    if sum(propProf)>0
        myFun(x) = LinearShoalBreak_Residual(x, H1[propProf], T1[propProf], DIR1[propProf], h1[propProf], ANGbati[propProf], Bcoef)
        # ub=2.0.*h1[2].*ones(size(h2l0[propProf]))
        # lb=zeros(size(h2l0[propProf]))
        # aux = nlboxsolve(myFun,h2l0[propProf],xtol=1e-3,method=:jfnk)
        # h2l=JFNK_ad(myFun,h2l0[propProf],lb,ub)
        
        # h2l=aux.zero
        # h2l=JFNK_new(myFun,h2l0[propProf],1e-5,200)
        h2l=optimize.newton_krylov(myFun,h2l0[propProf] ; method="minres")
        # h2l=linsolve(myFun,h2l0[propProf],tol=1e-5)
        H2l, DIR2l = LinearShoalBreak_ResidualVOL(h2l, H1[propProf],T1[propProf], DIR1[propProf], h1[propProf], ANGbati[propProf], Bcoef);                
        H2[propProf] = H2l
        DIR2[propProf] = DIR2l
        h2[propProf] = h2l
    end
        
    #     return H2; DIR2; h2
    
    return H2, DIR2, h2
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

function cartesianDir2nauticalDir(cDir)
    #     ###########################################################################    
    #     # Cartesian convention with 0 in East & positive counterclockwise TO
    #     # Nautical convention with 0 in North & positive clockwise. 
    #     ###########################################################################    
    
    nDir = 90.0 .- cDir
    nDir[nDir.<0] = 360.0 .+ nDir[nDir.<0]
    
    return nDir
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

function GroupCelerity(L,T,h)
    #     ###########################################################################    
    #     # CELERITY GROUP
    #     # L: wave lenght.
    #     # T: wave period.
    #     # h: depth of wave conditions.
    #     ###########################################################################       
    
    c = L./T
    k = 2 .*pi./L
    N = 1.0.+ 2.0 .*k.*h./sinh.(2.0 .*k.*h)
    Cg = c ./ 2.0 .* N
    
    return Cg
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

function interp_lon(x,lon,xq,varargin)

    # INTERP_LON interpolates a set of longitude angles [in deg]
    #
    # Usage: out = interp_lon(x,lon,xq)
    #
    # x & lon are vectors of length N.  function evalutes longitude 
    # (in deg -180..180) at points xq using unwrap & interp1()
    #
    # to specify interpolation method used in interp1; use
    # out = interp_lon(x,lon,xq,METHOD)
    #
    # Written by D.G. Long; 27 Nov 2017 
    
    ulon=unwrap(lon*pi/180)*180/pi
    if nargin>3
      out=interp1(x,ulon,xq,varargin[1])
    else()
      out=interp1(x,ulon,xq)
    end
    out=out%360
    out[out.>180]=out[out.>180]-360
    return out
end

function LinearShoal(H1, T1, DIR1, h1, h2, ANGbati)
    #     ###########################################################################    
    #     # Wave shoaling & refraction applying linear theory with parallel; rectilinear bathymetry.
    #     #    
    #     # INPUT:
    #     # H1:        initial wave height.
    #     # T1:        wave period.
    #     # DIR1:      initial wave direction. Nautical convention.
    #     # h1:        initial depth of wave conditions.
    #     # h2:        final depth of wave conditions.
    #     # ANGbati:   bathymetry angle; the normal of the shoreline. Cartesian convention
    #     #
    #     # OUTPUT:
    #     # H2:        wave height during breaking. Wave period is assumed invariant due to linear theory.
    #     # DIR2:      wave direction during breaking. Nautical convention.
    #     ###########################################################################
    # h2=real.(h2)
    relDir1 = rel_angle_cartesian(nauticalDir2cartesianDir(DIR1),ANGbati)
    L1, _ =RelDisp(h1,T1)
    L2, _ =RelDisp(h2,T1)
    CG1 = GroupCelerity(L1,T1,h1)
    CG2 = GroupCelerity(L2,T1,h2)
    relDir2 = Snell_Law(L1,L2,relDir1)
    KS = sqrt.(CG1./CG2)
    KR = sqrt.(cos.(relDir1.*pi./180.)./cos.(relDir2.*pi./180.))
    H2 = H1.*KS.*KR
    DIR2 = cartesianDir2nauticalDir(abs_angle_cartesian(relDir2,ANGbati))
    
    return H2, DIR2
end

function LinearShoalBreak_Residual(h2l, H1, T1, DIR1, h1, ANGbati, Bcoef)

    H2l, DIR2l = LinearShoal(H1, T1, DIR1, h1, h2l, ANGbati)
    H2comp = h2l.*Bcoef
    res = H2l-H2comp
    # res=res
    return res
end

function LinearShoalBreak_ResidualVOL(h2l, H1, T1, DIR1, h1, ANGbati, Bcoef)

    H2l, DIR2l = LinearShoal(H1, T1, DIR1, h1, h2l, ANGbati)
    H2comp = h2l.*Bcoef
    res = H2l-H2comp
    # res=res
    return H2l, DIR2l
end


function MILLER_DEAN_yeq(hb,tb,dirb,depthb,surge,tide,ammsla,aslr,sl,sc,dev,model,theta,AY0,qold,q,dc,dt,dx,D50,Hberm)
    #    #####################
    #    # Miller & Dean 2004 Modified josseaaa    
    #    # INPUT:
    #    #
    #    # wavec : wave climate at breaking during the new time step & 2 previous conditions
    #    #   hb  : wave height at breaking
    #    #   tb  : wave period at breaking
    #    #   dirb: wave direction at breaking
    #    #   depthb  : breaking depth
    #    #   surge: storm surge
    #    #   tide: astronomical tide
    #    #   mmsl: monthly mean sea level NOT the anomaly
    #    #   slr : mean sea level rise
    #    # sl    : source of sediment from the alongshore
    #    # sc    : source of sediment from the crosshore 
    #    # cnst  : model constants
    #    #   dc  : depth of closure
    #    #   D50 : mean sediment grain size()
    #    #   Hberm: berm height
    #    #   Hs12: wave height overcome 12hours in a year
    #    #   Ts12: representative period asociated to the wave height overcome 12 hours in a year
    #    # calp  : calibration parameters
    #    #   AY0 : Miller&Dean maximum progadation position without forcings
    #    # model : interaction with processes to turn on ["mmsl",'slr',"alongshore"]
    #    #
    #    # DEPENDENCIAS:
    #    # BruunRule; wast
    #    ######################################################################   
    #    # model : interaction with processes to turn on ["mmsl",'slr',"alongshore"]
    
        AYmmsla = 0
        AYslr = 0
        AYlst = 0
        AYcss = 0
        AYdune = 0
        dc = 0.5.*(dc[2:end,:].+dc[1:end-1,:]); ## trapezoidal rule for the closure depth in m+1/2
        ammsla = 0.5.*(ammsla[2:end].+ammsla[1:end-1]) ## trapezoidal rule for the mmsla in m+1/2
        aslr = 0.5.*(aslr[2:end].+aslr[1:end-1]) ## trapezoidal rule for the slr in m+1/2
    
    if model["mmsl"] == 1.0
        AYmmsla = -BruunRule(dc[:,end],D50,Hberm,ammsla)
    end
    
    if model["slr"]==1.0
        AYslr = -BruunRule(dc[:,end],D50,Hberm,aslr)
    end
    
    if model["alongshore"]==1.0
        AYlst = dt.*(-1.0./dc[:,end].*theta.*(q[2:end].-q[1:end-1] .+ sl[2:end,end].-sl[1:end-1,end])./dx + -1. ./dc[:,end-2].*
        (1-theta).*(qold[2:end].-qold[1:end-1] .+ sl[2:end,end-2].-sl[1:end-1,end-2])./dx)
    end
    
    if model["crosshore"]==1.0
        AYcss = dt.*(-1.0./dc[:,end].*theta.*0.5.*(sc[2:end,end].+sc[1:end-1,end])./dx .+ -1. ./dc[:,end-2].*(1-theta).*0.5.*(sc[2:end,end-2]+sc[1:end-1,end-2])./dx)
    end
    
    if model["dune"]==1.0
        AYdune = +1.0./dc[:,end].*dev
        AYdune[dev.<0] .= 0  # # guarantees when the storm finishes there is no progradation
    end     
                
        
    AY0 .= AY0 .+ AYmmsla .+ AYslr .+ AYlst .+ AYcss .+ AYdune
    wl = 0.106.*hb.+surge.+tide; # # total water level
    wl = 0.5.*(wl[2:end].+wl[1:end-1]); # # trapezoidal rule for the water level in m+1/2
    depthb = 0.5.*(depthb[2:end].+depthb[1:end-1]);  # # trapezoidal rule for the breaking depth in m+1/2    
    Wast = wast(depthb,D50)
    yeq = AY0 .- Wast.*(wl)./(Hberm.+depthb)
    
    return yeq, AY0    
end

function nauticalDir2cartesianDir(nDir)

    #     ###########################################################################    
    #     # Nautical convention with 0 in North & positive clockwise TO 
    #     # Cartesian convention with 0 in East & positive counterclockwise.
    #     ###########################################################################    
    # caso = 0
    # if isinstance[nDir,float] || isinstance[nDir,int]
    #     nDir = np.asarray[[nDir]]
    #     caso = 1
    # elseif isinstance[nDir,list]
    #     nDir = np.asarray[nDir]
    #     caso = 2
    # end
    # nDir=vec(nDir)
    cDir = 90.0 .- nDir
    for i in eachindex(cDir)
        if cDir[i] < -180.0
            cDir[i] = 360.0 .+ cDir[i]
        end
    end
    # if caso .== 1
    #     cDir = cDir.tolist[][0]
    # elseif caso .== 2
    #     cDir = cDir.tolist[]
    # end
            
    #     return cDir
    
    
    return cDir
end

function nauticalDir2cartesianDirP(nDir)

    #     ###########################################################################    
    #     # Nautical convention with 0 in North & positive clockwise TO 
    #     # Cartesian convention with 0 in East & positive counterclockwise.
    #     ###########################################################################    
    # caso = 0
    # if isinstance[nDir,float] || isinstance[nDir,int]
    #     nDir = np.asarray[[nDir]]
    #     caso = 1
    # elseif isinstance[nDir,list]
    #     nDir = np.asarray[nDir]
    #     caso = 2
    # end
    # nDir=vec(nDir)
    cDir = 90.0 - nDir
    
    if cDir < -180.0
        cDir = cDir + 360.0
    end
    
    # if caso .== 1
    #     cDir = cDir.tolist[][0]
    # elseif caso .== 2
    #     cDir = cDir.tolist[]
    # end
            
    #     return cDir
    
    
    return cDir
end

function pol2cart(rho,phi)

    x = rho .* cos(phi)
    y = rho .* sin(phi)
    return x,y
end

function rel_angle_cartesian(waveD, batiD)

    #     ###########################################################################    
    #     # Relative angle (in degrees) between wave direction & bathymetry with 
    #     # angles in cartesian coordinates, angle between [180,-180], 
    #     # 0 is in EAST & positive counterclockwise.
    #     #
    #     # INPUT:
    #     # waveD:    wave angle in Cartesian notation.
    #     # batiD:    bathymetry angle (normal to the shoreline) in Cartesian notation.
    #     #
    #     # OUTPUT:
    #     # relD:     relative wave angle between wave and bathymetry; 0 is the bathymetry & positive counterclockwise.
    #     ###########################################################################    
    # caso = 0
    # if isinstance[waveD,float] | isinstance[waveD,int]:
    #     waveD = np.asarray[[waveD]]
    #     caso = 1
    # elseif isinstance[waveD,list]
    #     waveD = np.asarray[waveD]
    #     caso = 2
    # end
    
    relD = waveD .- batiD
    relD[relD .> 180.] .= relD[relD .> 180.] .- 360.0
    relD[relD .< -180.] .= relD[relD .< -180.] .+ 360.0
    

    # if caso .== 1
    #     relD = relD.tolist[][0]
    # elseif caso .== 2
    #     relD = relD.tolist[]
    # end
       
    #     return relD
    return relD
end

function rel_angle_cartesianP(waveD, batiD)

    #     ###########################################################################    
    #     # Relative angle (in degrees) between wave direction & bathymetry with 
    #     # angles in cartesian coordinates, angle between [180,-180], 
    #     # 0 is in EAST & positive counterclockwise.
    #     #
    #     # INPUT:
    #     # waveD:    wave angle in Cartesian notation.
    #     # batiD:    bathymetry angle (normal to the shoreline) in Cartesian notation.
    #     #
    #     # OUTPUT:
    #     # relD:     relative wave angle between wave and bathymetry; 0 is the bathymetry & positive counterclockwise.
    #     ###########################################################################    
    # caso = 0
    # if isinstance[waveD,float] | isinstance[waveD,int]:
    #     waveD = np.asarray[[waveD]]
    #     caso = 1
    # elseif isinstance[waveD,list]
    #     waveD = np.asarray[waveD]
    #     caso = 2
    # end
    
    relD = waveD - batiD
    if relD > 180
        relD = relD - 360.0
    elseif relD < -180.0
        relD = relD + 360.0
    end
    

    # if caso .== 1
    #     relD = relD.tolist[][0]
    # elseif caso .== 2
    #     relD = relD.tolist[]
    # end
       
    #     return relD
    return relD
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

function residualLD_MD(ynew,y,r,dev,dt,dx,ti,hs,tp,dire,depth,hb,tb,dirb,depthb,surge,tide,mmsl,slr,q,
    yeq,req,Ts,D50,Hberm,DCT,AY0,kal,kacr,kero,trs,sc,sl,slope,dunezz,model,theta)

     # ynew[isnan.(ynew)].= 0.0
    # ynew[isinf.(ynew)].= 0.0
    X0 = trs["X0"]
    Y0 = trs["Y0"]
    phi = trs["phi"]
    # #    sl = bdc.ql
    # #    sc = bdc.qc
        
    crossp = 0.
    alongp = 0.
    mmslp = 0.
    slrp = 0.
    dunep = 0.
    
    
    if model["mmsl"] == 1.
        mmslp = 1.
    end
    if model["slr"] == 1.
        slrp=1.
    end
    if model["alongshore"] == 1.
        alongp=1.
    end
    if model["crosshore"] == 1.
        crossp=1.
    end
    if model["dune"] == 1.
        dunep=1.
    end 
    
              
    #     # calculating the shore positions
    XN, YN = abs_pos(X0,Y0,phi.*pi./180.,ynew)
    
    #     # calculating the shore angles
    alfas = zeros(size(hs,1),1)
    # alfam = np.zeros((hs.shape[1],),dtype=str)
    # alfas[1:-1],alfam[1:-1] = shore_angle[XN,YN,dire[:,ti]]
    alfas[2:end-1] = shore_angle(XN,YN,dire[:,ti])
    alfas[1] = alfas[2]; alfas[end] = alfas[end-1]; # # ghost condition for the relative angle()
    
    angulo_rel=-90; # # for Tairua is -90 according to how the reference system was defined
    #     # propagating waves with linear theory
    # #    hb[:,ti],dirb[:,ti],depthb[:,ti] = BreakingPropagation[hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas+90,"spectral"]#in Oregon this was the notation for getting the perpendicular to the shoreline
    hb[:,ti], dirb[:,ti], depthb[:,ti] = BreakingPropagation(hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas.-angulo_rel) #modificado LFP
    
    # #    hb[:,ti],dirb[:,ti],depthb[:,ti] = BreakingPropagationIntToRed[hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas+90,"spectral"]
    #     # calculating the alongshore sediment transport        
    #     sedtr = "cerq"
    
    kal = kal[:,ti]
    
    # #    kal = calp.K1[:,ti]
    # #    sedtr = "spm"
    # #    kal = calp.K[:,ti]
    
    dc = 0.5.*(DCT[2:end,ti-1:ti+1]+DCT[1:end-1,ti-1:ti+1]); # # trapezoidal rule for the closure depth in m+1/2
    
    # #    q[:,ti] = ALST[hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas+90,sedtr,kal]#in Oregon this was the notation for getting the perpendicular to the shoreline
    
    q[:,ti], _ = ALST(hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas.-angulo_rel,kal) # # for Tairua is -90 according to how the reference system was defined
    
    
    twl = dunep.*(RU2_Stockdon2006(slope[:,ti],hs[:,ti],tp[:,ti])+0.5.*(surge[2:end,ti]+surge[1:end-1,ti])+0.5.*
    (tide[2:end,ti]+tide[1:end-1,ti])+mmslp*(0.5.*(mmsl[2:end,ti]+mmsl[1:end-1,ti]))+slrp.*(0.5.*(slr[2:end,ti]+slr[1:end-1,ti]))); ##TWL RUGGIERO
    
    
    twl[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1]]=dunezz[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1],2]
    # # we just run collapsing events; not overwash

    #     # calculating the crossshore equilibrium position

    yeq[:,ti], AY0[:,ti] = MILLER_DEAN_yeq(hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],surge[:,ti],tide[:,ti],mmslp.*
    (mmsl[:,ti].-mmsl[:,ti-1]),slrp.*(slr[:,ti].-slr[:,ti-1]),sl[:,ti-1:ti+1],sc[:,ti-1:ti+1],dunep.*(dev[:,ti]),model,
    theta,AY0[:,ti-1],q[:,ti-1],q[:,ti],DCT[:,ti-1:ti+1],dt,dx,D50[:,ti],Hberm[:,ti])
        
    #     # defining if the transect is eroding | accreting due to cross-shore forcings.
    rcross = theta.*(yeq[:,ti].-ynew)+(1-theta).*(yeq[:,ti-1].-y)
    kcr = zeros(size(rcross))
    kcr[rcross .>= 0] .= kacr[vec(rcross .>= 0),ti]
    kcr[rcross .< 0] .= kero[vec(rcross) .< 0,ti]
    # # lets do it proportional to Hb2
    # #    kcr[(hb[ti,1:]).>0.05]=kcr[(hb[ti,1:]).>0.05]*(0.5*(hb[ti,1:]+hb[ti,:-1])**2)[(hb[ti,1:]).>0.05]
        
    # kcr=kcr.*(0.5.*(hb[2:end,ti]+hb[1:end-1,ti]).^3)

    # q[1,:].=q[2,:]
    # q[end,:].=q[end-1,:]
    # kcr[1]=kcr[2]
    # kcr[end]=kcr[end-1]
    # yeq[1,:].=kcr[2,:]
    # yeq[end,:].=yeq[end-1,:]


    return (ynew.-y)./dt .+ alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]) .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(-kcr.*(theta.*(yeq[:,ti].-ynew).+(1-theta).*(yeq[:,ti-1].-y))) .+ theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1]

end

function residualLDVOL_MD(ynew,y,r,dev,dt,dx,ti,hs,tp,dire,depth,hb,tb,dirb,depthb,surge,tide,mmsl,slr,q,
    yeq,req,Ts,D50,Hberm,DCT,AY0,kal,kacr,kero,trs,sc,sl,slope,dunezz,model,theta)
        
    X0 = trs["X0"]
    Y0 = trs["Y0"]
    phi = trs["phi"]
    # #    sl = bdc.ql
    # #    sc = bdc.qc
        
    crossp = 0.
    alongp = 0.
    mmslp = 0.
    slrp = 0.
    dunep = 0.
    
    
    if model["mmsl"] == 1.
        mmslp = 1.
    end
    if model["slr"] == 1.
        slrp=1.
    end
    if model["alongshore"] == 1.
        alongp=1.
    end
    if model["crosshore"] == 1.
        crossp=1.
    end
    if model["dune"] == 1.
        dunep=1.
    end 
    
              
    #     # calculating the shore positions
    XN, YN = abs_pos(X0,Y0,phi.*pi./180.,ynew)
    
    #     # calculating the shore angles
    alfas = zeros(size(hs,1),1)
    # alfam = np.zeros((hs.shape[1],),dtype=str)
    # alfas[1:-1],alfam[1:-1] = shore_angle[XN,YN,dire[:,ti]]
    alfas[2:end-1] = shore_angle(XN,YN,dire[:,ti])
    alfas[1] = alfas[2]; alfas[end] = alfas[end-1]; # # ghost condition for the relative angle()
    
    angulo_rel=-90; # # for Tairua is -90 according to how the reference system was defined
    #     # propagating waves with linear theory
    # #    hb[:,ti],dirb[:,ti],depthb[:,ti] = BreakingPropagation[hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas+90,"spectral"]#in Oregon this was the notation for getting the perpendicular to the shoreline
    hb[:,ti], dirb[:,ti], depthb[:,ti] = BreakingPropagation(hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas.-angulo_rel) #modificado LFP
    
    #     sedtr = "cerq"
    
    kal = kal[:,ti]
    
        
    dc = 0.5.*(DCT[2:end,ti-1:ti+1]+DCT[1:end-1,ti-1:ti+1]) # # trapezoidal rule for the closure depth in m+1/2
    
    # #    q[:,ti] = ALST[hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas+90,sedtr,kal]#in Oregon this was the notation for getting the perpendicular to the shoreline
    
    q[:,ti], q0 = ALST(hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas.-angulo_rel,kal) # # for Tairua is -90 according to how the reference system was defined
    
    
    twl = dunep.*(RU2_Stockdon2006(slope[:,ti],hs[:,ti],tp[:,ti])+0.5.*(surge[2:end,ti]+surge[1:end-1,ti])+0.5.*
    (tide[2:end,ti]+tide[1:end-1,ti])+mmslp*(0.5.*(mmsl[2:end,ti]+mmsl[1:end-1,ti]))+slrp.*(0.5.*(slr[2:end,ti]+slr[1:end-1,ti]))); ##TWL RUGGIERO
    
    
    twl[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1]]=dunezz[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1],2]
    
    #     # calculating the crossshore equilibrium position

    yeq[:,ti], AY0[:,ti] = MILLER_DEAN_yeq(hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],surge[:,ti],tide[:,ti],mmslp.*
    (mmsl[:,ti].-mmsl[:,ti-1]),slrp.*(slr[:,ti].-slr[:,ti-1]),sl[:,ti-1:ti+1],sc[:,ti-1:ti+1],dunep.*(dev[:,ti]),model,
    theta,AY0[:,ti-1],q[:,ti-1],q[:,ti],DCT[:,ti-1:ti+1],dt,dx,D50[:,ti],Hberm[:,ti])
        
    #     # defining if the transect is eroding | accreting due to cross-shore forcings.
    rcross = theta.*(yeq[:,ti].-ynew)+(1-theta).*(yeq[:,ti-1].-y)
    kcr = zeros(size(rcross))
    kcr[rcross .>= 0] .= kacr[vec(rcross .>= 0),ti]
    kcr[rcross .< 0] .= kero[vec(rcross) .< 0,ti]
        # # lets do it proportional to Hb2
        # #    kcr[(hb[ti,1:]).>0.05]=kcr[(hb[ti,1:]).>0.05]*(0.5*(hb[ti,1:]+hb[ti,:-1])**2)[(hb[ti,1:]).>0.05]
            
    # kcr=kcr.*(0.5.*(hb[2:end,ti]+hb[1:end-1,ti]).^3)
    
    # q[1,:].=q[2,:]
    # q[end,:].=q[end-1,:]
    # kcr[1] = kcr[2]
    # kcr[end] = kcr[end-1]
    # yeq[1,:].=kcr[2,:]
    # yeq[end,:].=yeq[end-1,:]
    # # lets do it proportional to Hb2
    # #    kcr[(hb[ti,1:]).>0.05]=kcr[(hb[ti,1:]).>0.05]*(0.5*(hb[ti,1:]+hb[ti,:-1])**2)[(hb[ti,1:]).>0.05]
        
    # # sediment budget Northerly transport positive
    sedbgt = -1.0*(alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1])
     .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(-kcr.*
     (theta.*(yeq[:,ti].-ynew).+(1-theta).*(yeq[:,ti-1].-y))) + theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*
     (sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1])
                
    sedbgtal = -1.0*(alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end].+(1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]))
                
    sedbgtcr = -1.0*(crossp.*(-kcr.*(theta.*(yeq[:,ti].-ynew)+(1-theta).*(yeq[:,ti-1].-y))))
    
    sedbgtbc = -1.0*(theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ 
    theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1])
    
    resid = (ynew.-y)./dt + alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]) 
        .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(-kcr.*(theta.*
        (yeq[:,ti].-ynew)+(1-theta).*(yeq[:,ti-1].-y))) + theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1]
    
    
    
    ctAL = dt.*(-alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]))
    ctCS = dt.*(crossp.*(kcr.*(theta.*(yeq[:,ti].-ynew)+(1-theta).*(yeq[:,ti-1].-y))))
    
    # Ylt = ctAL./(abs.(ctAL).+abs.(ctCS)) .* ynew
    # Yct = ctCS./(abs.(ctAL).+abs.(ctCS)) .* ynew

    kcerc = 0.5.*(kal[2:end]+kal[1:end-1])

    dQdx = (theta.*(q[2:end,ti].-q[1:end-1,ti])./dx .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx)./kcerc
    # dQdx[1]=0
    # dQdx[end]=0
    # dQdx = (q[2:end,ti].-q[1:end-1,ti])./dx
    # resid=resid'
    
    # rti=r[:,ti]
    # devti=dev[:,ti]
    # hbti=hb[:,ti]
    # dirbti=dirb[:,ti]
    # depthbti=depthb[:,ti]
    # qti=q[:,ti]
    # yeqti=yeq[:,ti]
    # AY0ti=AY0[:,ti]
    # reqti=req[:,ti]
    # Tsti=Ts[:,ti]
    
    #    return resid, rti, devti, hbti,dirbti, depthbti, twl, qti, yeqti, AY0ti,reqti, Tsti, sedbgt, sedbgtal,sedbgtcr, sedbgtbc 
    return resid, r[:,ti], dev[:,ti], hb[:,ti], dirb[:,ti], depthb[:,ti], twl, q[:,ti],q0, yeq[:,ti], AY0[:,ti], req[:,ti], Ts[:,ti], sedbgt, sedbgtal, sedbgtcr, sedbgtbc, dQdx, ctAL,ctCS
    # #    return resid, XN, YN, alfas, alfam, hb[-1,:], dirb[-1,:], depthb[-1,:], q[-1,:], yeq[-1,:], AY0[-1,:]

end

function residualLD_Y09(ynew,y,r,dev,dt,dx,ti,hs,tp,dire,depth,hb,tb,dirb,depthb,surge,tide,mmsl,slr,q,
    yeq,req,Ts,D50,Hberm,DCT,DY,kal,HSB,DT,trs,sc,sl,slope,dunezz,model,theta, tau,r_Brunn)

    ynew=real.(ynew)
    # ynew[isnan.(ynew)].= 0.0
    # ynew[isinf.(ynew)].= 0.0
    X0 = trs["X0"]
    Y0 = trs["Y0"]
    phi = trs["phi"]
    # #    sl = bdc.ql
    # #    sc = bdc.qc
        
    crossp = 0.
    alongp = 0.
    mmslp = 0.
    slrp = 0.
    dunep = 0.
    
    
    if model["mmsl"] == 1.
        mmslp = 1.
    end
    if model["slr"] == 1.
        slrp=1.
    end
    if model["alongshore"] == 1.
        alongp=1.
    end
    if model["crosshore"] == 1.
        crossp=1.
    end
    if model["dune"] == 1.
        dunep=1.
    end 
    
              
    #     # calculating the shore positions
    XN, YN = abs_pos(X0,Y0,phi.*pi./180.,ynew)
    
    #     # calculating the shore angles
    alfas = zeros(size(hs,1),1)
    # alfam = np.zeros((hs.shape[1],),dtype=str)
    # alfas[1:-1],alfam[1:-1] = shore_angle[XN,YN,dire[:,ti]]
    alfas[2:end-1] = shore_angle(XN,YN,dire[:,ti])
    alfas[1] = alfas[2]; alfas[end] = alfas[end-1]; # # ghost condition for the relative angle()
    
    angulo_rel=-90; # # for Tairua is -90 according to how the reference system was defined
    #     # propagating waves with linear theory
    # #    hb[:,ti],dirb[:,ti],depthb[:,ti] = BreakingPropagation[hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas+90,"spectral"]#in Oregon this was the notation for getting the perpendicular to the shoreline
    hb[:,ti], dirb[:,ti], depthb[:,ti] = BreakingPropagation(hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas.-angulo_rel) #modificado LFP
    
    # #    hb[:,ti],dirb[:,ti],depthb[:,ti] = BreakingPropagationIntToRed[hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas+90,"spectral"]
    #     # calculating the alongshore sediment transport        
    #     sedtr = "cerq"
    
    kal = kal[:,ti]
    
    # #    kal = calp.K1[:,ti]
    # #    sedtr = "spm"
    # #    kal = calp.K[:,ti]
    
    dc = 0.5.*(DCT[2:end,ti-1:ti+1]+DCT[1:end-1,ti-1:ti+1]); # # trapezoidal rule for the closure depth in m+1/2
    
    # #    q[:,ti] = ALST[hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas+90,sedtr,kal]#in Oregon this was the notation for getting the perpendicular to the shoreline
    
    q[:,ti], _ = ALST(hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas.-angulo_rel,kal) # # for Tairua is -90 according to how the reference system was defined
    
    
    twl = dunep.*(RU2_Stockdon2006(slope[:,ti],hs[:,ti],tp[:,ti])+0.5.*(surge[2:end,ti]+surge[1:end-1,ti])+0.5.*
    (tide[2:end,ti]+tide[1:end-1,ti])+mmslp*(0.5.*(mmsl[2:end,ti]+mmsl[1:end-1,ti]))+slrp.*(0.5.*(slr[2:end,ti]+slr[1:end-1,ti]))); ##TWL RUGGIERO
    
    
    twl[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1]]=dunezz[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1],2]
    # # we just run collapsing events; not overwash

    #     # calculating the crossshore equilibrium position
    hss = 0.5.*(hs[2:end,ti]+hs[1:end-1,ti])
    hss[1] = hss[2]
    hss[end] = hss[end-1] 
    yeq[:,ti], tau[:,ti] = Y09_MOD(hss,HSB[:,ti],DY[:,ti],DT[:,ti])
    kcr = 1 ./tau[:,ti]

    slrt=0.5.*(slr[2:end,:].+slr[1:end-1,:])
    r_Brunn[:,ti] = BruunRule(dc[:,ti],D50[:,ti],Hberm[:,ti],slrt[:,ti].-slrt[:,ti-1])
    # Azin=ADEAN(D50)
    # tbeta=dc[5,1]*(0.55*Azin[5,1])^1.5/hb[5,ti]
    # q[1,:].=q[2,:]
    # q[end,:].=q[end-1,:]
    # kcr[1]=kcr[2]
    # kcr[end]=kcr[end-1]
    # yeq[1,:].=kcr[2,:]
    # yeq[end,:].=yeq[end-1,:]

    return (ynew.-y)./dt .+ alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]) .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(-kcr.*(theta.*(yeq[:,ti].-ynew).+(1-theta).*(yeq[:,ti-1].-y))) .+ theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1] .+ (theta.*r_Brunn[:,ti] .+ (1-theta).*r_Brunn[:,ti-1])./dt

end

function residualLDVOL_Y09(ynew,y,r,dev,dt,dx,ti,hs,tp,dire,depth,hb,tb,dirb,depthb,surge,tide,mmsl,slr,q,
    yeq,req,Ts,D50,Hberm,DCT,DY,kal,HSB,DT,trs,sc,sl,slope,dunezz,model,theta,tau,r_Brunn)
        
    X0 = trs["X0"]
    Y0 = trs["Y0"]
    phi = trs["phi"]
    # #    sl = bdc.ql
    # #    sc = bdc.qc
        
    crossp = 0.
    alongp = 0.
    mmslp = 0.
    slrp = 0.
    dunep = 0.
    
    
    if model["mmsl"] == 1.
        mmslp = 1.
    end
    if model["slr"] == 1.
        slrp=1.
    end
    if model["alongshore"] == 1.
        alongp=1.
    end
    if model["crosshore"] == 1.
        crossp=1.
    end
    if model["dune"] == 1.
        dunep=1.
    end 
    
              
    #     # calculating the shore positions
    XN, YN = abs_pos(X0,Y0,phi.*pi./180.,ynew)
    
    #     # calculating the shore angles
    alfas = zeros(size(hs,1),1)
    # alfam = np.zeros((hs.shape[1],),dtype=str)
    # alfas[1:-1],alfam[1:-1] = shore_angle[XN,YN,dire[:,ti]]
    alfas[2:end-1] = shore_angle(XN,YN,dire[:,ti])
    alfas[1] = alfas[2]; alfas[end] = alfas[end-1]; # # ghost condition for the relative angle()
    
    angulo_rel=-90; # # for Tairua is -90 according to how the reference system was defined
    #     # propagating waves with linear theory
    # #    hb[:,ti],dirb[:,ti],depthb[:,ti] = BreakingPropagation[hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas+90,"spectral"]#in Oregon this was the notation for getting the perpendicular to the shoreline
    hb[:,ti], dirb[:,ti], depthb[:,ti] = BreakingPropagation(hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas.-angulo_rel) #modificado LFP
    
    #     sedtr = "cerq"
    
    kal = kal[:,ti]
    
        
    dc = 0.5.*(DCT[2:end,ti-1:ti+1]+DCT[1:end-1,ti-1:ti+1]) # # trapezoidal rule for the closure depth in m+1/2
    
    # #    q[:,ti] = ALST[hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas+90,sedtr,kal]#in Oregon this was the notation for getting the perpendicular to the shoreline
    
    q[:,ti], q0 = ALST(hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas.-angulo_rel,kal) # # for Tairua is -90 according to how the reference system was defined
    
    
    twl = dunep.*(RU2_Stockdon2006(slope[:,ti],hs[:,ti],tp[:,ti])+0.5.*(surge[2:end,ti]+surge[1:end-1,ti])+0.5.*
    (tide[2:end,ti]+tide[1:end-1,ti])+mmslp*(0.5.*(mmsl[2:end,ti]+mmsl[1:end-1,ti]))+slrp.*(0.5.*(slr[2:end,ti]+slr[1:end-1,ti]))); ##TWL RUGGIERO
    
    
    twl[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1]]=dunezz[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1],2]
    
    #     # calculating the crossshore equilibrium position
    
    hss = 0.5.*(hs[2:end,ti]+hs[1:end-1,ti])
    # hss[1] = hss[2]
    # hss[end] = hss[end-1] 
    yeq[:,ti], tau[:,ti] = Y09_MOD(hss,HSB[:,ti],DY[:,ti],DT[:,ti])
    kcr = 1 ./tau[:,ti]
    slrt=0.5.*(slr[2:end,:].+slr[1:end-1,:])
    r_Brunn[:,ti] = BruunRule(dc[:,ti],D50[:,ti],Hberm[:,ti],slrt[:,ti].-slrt[:,ti-1])
    # Azin=ADEAN(D50)
    # tbeta=dc[5,1]*(0.55*Azin[5,1])^1.5/hb[5,ti]
    # q[1,:].=q[2,:]
    # q[end,:].=q[end-1,:]
    # kcr[1]=kcr[2]
    # kcr[end]=kcr[end-1]
    # yeq[1,:].=yeq[2,:]
    # yeq[end,:].=yeq[end-1,:]
    # # lets do it proportional to Hb2
    # #    kcr[(hb[ti,1:]).>0.05]=kcr[(hb[ti,1:]).>0.05]*(0.5*(hb[ti,1:]+hb[ti,:-1])**2)[(hb[ti,1:]).>0.05]
        
    # # sediment budget Northerly transport positive
    sedbgt = -1.0*(alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1])
     .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(-kcr.*
     (theta.*(yeq[:,ti].-ynew).+(1-theta).*(yeq[:,ti-1].-y))) + theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*
     (sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1])
                
    sedbgtal = -1.0*(alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end].+(1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]))
                
    sedbgtcr = -1.0*(crossp.*(-kcr.*(theta.*(yeq[:,ti].-ynew)+(1-theta).*(yeq[:,ti-1].-y))))
    
    sedbgtbc = -1.0*(theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ 
    theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1])
    
    resid = (ynew.-y)./dt + alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]) 
        .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(-kcr.*(theta.*
        (yeq[:,ti].-ynew)+(1-theta).*(yeq[:,ti-1].-y))) + theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1].+
        (theta.*r_Brunn[:,ti] .+ (1-theta).*r_Brunn[:,ti-1])./50.0./dt
    
    
    
    ctAL = dt.*(-alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]))
    ctCS = dt.*(kcr.*(theta.*(yeq[:,ti].-ynew)+(1-theta).*(yeq[:,ti-1].-y)))
    
    # Ylt = ctAL./(abs.(ctAL).+abs.(ctCS)) .* ynew
    # Yct = ctCS./(abs.(ctAL).+abs.(ctCS)) .* ynew

    kcerc = 0.5.*(kal[2:end]+kal[1:end-1])

    dQdx =(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx)./kcerc

    
    
    #    return resid, rti, devti, hbti,dirbti, depthbti, twl, qti, yeqti, AY0ti,reqti, Tsti, sedbgt, sedbgtal,sedbgtcr, sedbgtbc 
    return resid, r[:,ti], dev[:,ti], hb[:,ti], dirb[:,ti], depthb[:,ti], twl, q[:,ti],q0, yeq[:,ti], req[:,ti], Ts[:,ti], sedbgt, sedbgtal, sedbgtcr, sedbgtbc, dQdx, ctAL,ctCS, tau[:,ti], r_Brunn[:,ti]
    # #    return resid, XN, YN, alfas, alfam, hb[-1,:], dirb[-1,:], depthb[-1,:], q[-1,:], yeq[-1,:], AY0[-1,:]

end

function residualLD_SF(ynew,y,r,dev,dt,dx,ti,hs,tp,dire,depth,hb,tb,dirb,depthb,surge,tide,mmsl,slr,q,
    req,Ts,D50,Hberm,DCT,kal,trs,sc,sl,slope,dunezz,model,theta,F,OmegaEQ,cacr_SF,cero_SF,Omega,dpO, r_Brunn)

    ynew=real.(ynew)
    # ynew[isnan.(ynew)].= 0.0
    # ynew[isinf.(ynew)].= 0.0
    X0 = trs["X0"]
    Y0 = trs["Y0"]
    phi = trs["phi"]
    # #    sl = bdc.ql
    # #    sc = bdc.qc
        
    crossp = 0.
    alongp = 0.
    mmslp = 0.
    slrp = 0.
    dunep = 0.
    
    
    if model["mmsl"] == 1.
        mmslp = 1.
    end
    if model["slr"] == 1.
        slrp=1.
    end
    if model["alongshore"] == 1.
        alongp=1.
    end
    if model["crosshore"] == 1.
        crossp=1.
    end
    if model["dune"] == 1.
        dunep=1.
    end 
    
              
    #     # calculating the shore positions
    XN, YN = abs_pos(X0,Y0,phi.*pi./180.,ynew)
    
    #     # calculating the shore angles
    alfas = zeros(size(hs,1),1)
    # alfam = np.zeros((hs.shape[1],),dtype=str)
    # alfas[1:-1],alfam[1:-1] = shore_angle[XN,YN,dire[:,ti]]
    alfas[2:end-1] = shore_angle(XN,YN,dire[:,ti])
    alfas[1] = alfas[2]; alfas[end] = alfas[end-1]; # # ghost condition for the relative angle()
    
    angulo_rel=-90; # # for Tairua is -90 according to how the reference system was defined
    #     # propagating waves with linear theory
    # #    hb[:,ti],dirb[:,ti],depthb[:,ti] = BreakingPropagation[hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas+90,"spectral"]#in Oregon this was the notation for getting the perpendicular to the shoreline
    hb[:,ti], dirb[:,ti], depthb[:,ti] = BreakingPropagation(hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas.-angulo_rel) #modificado LFP
    
    
    kal = kal[:,ti]
    
    # #    kal = calp.K1[:,ti]
    # #    sedtr = "spm"
    # #    kal = calp.K[:,ti]
    
    dc = 0.5.*(DCT[2:end,ti-1:ti+1]+DCT[1:end-1,ti-1:ti+1]); # # trapezoidal rule for the closure depth in m+1/2
    
    # #    q[:,ti] = ALST[hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas+90,sedtr,kal]#in Oregon this was the notation for getting the perpendicular to the shoreline
    
    q[:,ti], _ = ALST(hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas.-angulo_rel,kal) # # for Tairua is -90 according to how the reference system was defined
    
    
    twl = dunep.*(RU2_Stockdon2006(slope[:,ti],hs[:,ti],tp[:,ti])+0.5.*(surge[2:end,ti]+surge[1:end-1,ti])+0.5.*
    (tide[2:end,ti]+tide[1:end-1,ti])+mmslp*(0.5.*(mmsl[2:end,ti]+mmsl[1:end-1,ti]))+slrp.*(0.5.*(slr[2:end,ti]+slr[1:end-1,ti]))); ##TWL RUGGIERO
    
    
    twl[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1]]=dunezz[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1],2]
    # # we just run collapsing events; not overwash

    #     # calculating the crossshore equilibrium position

    F[:,ti], _ = ShoreFor(OmegaEQ,dpO,tb[:,ti],hb[:,ti],depthb[:,ti],D50[:,ti],Omega,dt)
    F_c=copy(F)
    c_SF = zeros(size(cacr_SF))
    c_SF[vec(F[:,ti] .>= 0)] .= cacr_SF[vec(F[:,ti] .>= 0)]
    c_SF[vec(F[:,ti] .< 0)] .= cero_SF[vec(F[:,ti] .< 0)]
    # # lets do it proportional to Hb2
    # #    kcr[(hb[ti,1:]).>0.05]=kcr[(hb[ti,1:]).>0.05]*(0.5*(hb[ti,1:]+hb[ti,:-1])**2)[(hb[ti,1:]).>0.05]

    # for i in eachindex(F_c[1,:])
    #     posero =F_c[:,i].<0
    #     F_c[posero,i] .= F_c[posero,i].*r_SF[posero]
    # end

    slrt=0.5.*(slr[2:end,:].+slr[1:end-1,:])
    r_Brunn[:,ti] = BruunRule(dc[:,ti],D50[:,ti],Hberm[:,ti],slrt[:,ti])#.-slrt[:,ti-1]
    # Azin=ADEAN(D50)
    # tbeta=dc[5,1]*(0.55*Azin[5,1])^1.5/hb[5,ti]
    # q[1,:].=q[2,:]
    # q[end,:].=q[end-1,:]
    # c_SF[1]=c_SF[2]
    # c_SF[end]=c_SF[end-1]
    # F_c[1,:].=F_c[2,:]
    # F_c[end,:].=F_c[end-1,:]

    return (ynew.-y)./dt .+ alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]) .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(-c_SF.*(theta.*(F_c[:,ti]).+(1-theta).*(F_c[:,ti-1]))) .+ theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1]#.+ (theta.*r_Brunn[:,ti] .+ (1-theta).*r_Brunn[:,ti-1])./dt

end

function residualLDVOL_SF(ynew,y,r,dev,dt,dx,ti,hs,tp,dire,depth,hb,tb,dirb,depthb,surge,tide,mmsl,slr,q,
    req,Ts,D50,Hberm,DCT,kal,trs,sc,sl,slope,dunezz,model,theta,F,OmegaEQ,cacr_SF,cero_SF,Omega,dpO, r_Brunn)
        
    X0 = trs["X0"]
    Y0 = trs["Y0"]
    phi = trs["phi"]
    # #    sl = bdc.ql
    # #    sc = bdc.qc
        
    crossp = 0.
    alongp = 0.
    mmslp = 0.
    slrp = 0.
    dunep = 0.
    
    
    if model["mmsl"] == 1.
        mmslp = 1.
    end
    if model["slr"] == 1.
        slrp=1.
    end
    if model["alongshore"] == 1.
        alongp=1.
    end
    if model["crosshore"] == 1.
        crossp=1.
    end
    if model["dune"] == 1.
        dunep=1.
    end 
    
              
    #     # calculating the shore positions
    XN, YN = abs_pos(X0,Y0,phi.*pi./180.,ynew)
    
    #     # calculating the shore angles
    alfas = zeros(size(hs,1),1)
    # alfam = np.zeros((hs.shape[1],),dtype=str)
    # alfas[1:-1],alfam[1:-1] = shore_angle[XN,YN,dire[:,ti]]
    alfas[2:end-1] = shore_angle(XN,YN,dire[:,ti])
    alfas[1] = alfas[2]; alfas[end] = alfas[end-1]; # # ghost condition for the relative angle()
    
    angulo_rel=-90; # # for Tairua is -90 according to how the reference system was defined
    #     # propagating waves with linear theory
    # #    hb[:,ti],dirb[:,ti],depthb[:,ti] = BreakingPropagation[hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas+90,"spectral"]#in Oregon this was the notation for getting the perpendicular to the shoreline
    hb[:,ti], dirb[:,ti], depthb[:,ti] = BreakingPropagation(hs[:,ti],tp[:,ti],dire[:,ti],depth[:,ti],alfas.-angulo_rel) #modificado LFP
    
    #     sedtr = "cerq"
    
    kal = kal[:,ti]
    
        
    dc = 0.5.*(DCT[2:end,ti-1:ti+1]+DCT[1:end-1,ti-1:ti+1]) # # trapezoidal rule for the closure depth in m+1/2
    
    # #    q[:,ti] = ALST[hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas+90,sedtr,kal]#in Oregon this was the notation for getting the perpendicular to the shoreline
    
    q[:,ti], q0 = ALST(hb[:,ti],tb[:,ti],dirb[:,ti],depthb[:,ti],alfas.-angulo_rel,kal) # # for Tairua is -90 according to how the reference system was defined
    
    
    twl = dunep.*(RU2_Stockdon2006(slope[:,ti],hs[:,ti],tp[:,ti])+0.5.*(surge[2:end,ti]+surge[1:end-1,ti])+0.5.*
    (tide[2:end,ti]+tide[1:end-1,ti])+mmslp*(0.5.*(mmsl[2:end,ti]+mmsl[1:end-1,ti]))+slrp.*(0.5.*(slr[2:end,ti]+slr[1:end-1,ti]))); ##TWL RUGGIERO
    
    
    twl[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1]]=dunezz[twl-dunezz[:,1].>dunezz[:,2]-dunezz[:,1],2]
    
    #     # calculating the crossshore equilibrium position


    F[:,ti], Omega[:,ti] = ShoreFor(OmegaEQ,dpO,tb[:,ti],hb[:,ti],depthb[:,ti],D50[:,ti],Omega,dt)
    
    F_c=copy(F)
    c_SF = zeros(size(cacr_SF))
    c_SF[vec(F[:,ti] .>= 0)] .= cacr_SF[vec(F[:,ti] .>= 0)]
    c_SF[vec(F[:,ti] .< 0)] .= cero_SF[vec(F[:,ti] .< 0)]

    # F_c = copy(F)
    # for i in eachindex(F_c[1,:])
    #     posero =F_c[:,i].<0
    #     F_c[posero,i] .= F_c[posero,i].*r_SF[posero]
    # end

    slrt=0.5.*(slr[2:end,:].+slr[1:end-1,:])
    r_Brunn[:,ti] = BruunRule(dc[:,ti],D50[:,ti],Hberm[:,ti],slrt[:,ti].-slrt[:,ti-1])
    # Azin=ADEAN(D50)
    # tbeta=dc[5,1]*(0.55*Azin[5,1])^1.5/hb[5,ti]
    # q[1,:].=q[2,:]
    # q[end,:].=q[end-1,:]
    # c_SF[1]=c_SF[2]
    # c_SF[end]=c_SF[end-1]
    # F_c[1,:].=F_c[2,:]
    # F_c[end,:].=F_c[end-1,:]

    # F[1,:].=F[2,:]
    # F[end,:].=F[end-1,:]
    # Omega[1,:].=Omega[2,:]
    # Omega[end,:].=Omega[end-1,:]
    # # lets do it proportional to Hb2
    # #    kcr[(hb[ti,1:]).>0.05]=kcr[(hb[ti,1:]).>0.05]*(0.5*(hb[ti,1:]+hb[ti,:-1])**2)[(hb[ti,1:]).>0.05]
        
    # # sediment budget Northerly transport positive
    sedbgt = -1.0*(alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1])
     .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(c_SF.*
     (theta.*(F_c[:,ti]).+(1-theta).*(F_c[:,ti-1]))) + theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*
     (sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1])
                
    sedbgtal = -1.0*(alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end].+(1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]))
                
    sedbgtcr = -1.0*(crossp.*(-c_SF.*(theta.*(F_c[:,ti])+(1-theta).*(F_c[:,ti-1]))))
    
    sedbgtbc = -1.0*(theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ 
    theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1])
    
    resid = (ynew.-y)./dt + alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]) 
        .+ theta.*(sl[2:end,ti].-sl[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(sl[2:end,ti-1].-sl[1:end-1,ti-1])./dx./dc[:,end-1] .+ crossp.*(-c_SF.*(theta.*
        (F_c[:,ti])+(1-theta).*(F_c[:,ti-1]))) + theta.*0.5.*(sc[2:end,ti].+sc[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*0.5.*(sc[2:end,ti-1].+sc[1:end-1,ti-1])./dx./dc[:,end-1] .+
        (theta.*r_Brunn[:,ti] .+ (1-theta).*r_Brunn[:,ti-1])./50.0./dt
    
    
    
    ctAL = dt.*(-alongp.*(theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,end] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,end-1]))
    ctCS = dt.*(crossp.*(c_SF.*(theta.*(F_c[:,ti])+(1-theta).*(F_c[:,ti-1]))))
    
    # Ylt = ctAL./(abs.(ctAL).+abs.(ctCS)) .* ynew
    # Yct = ctCS./(abs.(ctAL).+abs.(ctCS)) .* ynew

    kcerc = 0.5.*(kal[2:end]+kal[1:end-1])

    dQdx = (theta.*(q[2:end,ti].-q[1:end-1,ti])./dx .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx)./kcerc

    
    return resid, r[:,ti], dev[:,ti], hb[:,ti], dirb[:,ti], depthb[:,ti], twl, q[:,ti],q0, req[:,ti], Ts[:,ti], sedbgt, sedbgtal, sedbgtcr, sedbgtbc, dQdx, ctAL,ctCS, F[:,ti], Omega[:,ti], OmegaEQ, r_Brunn[:,ti]

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

function shore_angle(XN,YN,wave_angle)
    #     #####################    
    #     # INPUT:
    #     #
    #     # XN : x coordinate
    #     # YN : y coordinate
    #     # wave_angle : wave angle()
    #     #
    #     # OUTPUT:
    #     # shoreAng : angle of the shoreline to compute sediment transport
    #     #####################
    
    shoreAng = zeros(1,length(XN)-1)
    # method = zeros(length(XN)-1)
    
    for i in eachindex(shoreAng)
        alfa = atan(YN[i+1]-YN[i],XN[i+1]-XN[i]).*180.0./pi
        relang = rel_angle_cartesianP(nauticalDir2cartesianDirP(wave_angle[i+1]),alfa)
    #         # check the relative angle between the waves & the shoreline orientation [not the normal to the shoreline]
        if abs(relang).>=45 && abs(relang).<=180-45
            shoreAng[i] = alfa
    #         method[i] = "CenDiff"
        elseif abs(relang).<45  && i.<length(shoreAng)-1
            try
                shoreAng[i] = atan(YN[i+2]-YN[i+1],XN[i+2]-XN[i+1]).*180.0./pi
    #             method[i] = "UpWind"
            catch
                shoreAng[i] = alfa
    #             method[i] = "CenDiff"
            end
                    
        elseif abs(relang).>180 - 45 && i.>1
            try
                shoreAng[i] = atan(YN[i]-YN[i-1],XN[i]-XN[i-1]).*180.0./pi
    #             method[i] = "UpWind"
            catch
                shoreAng[i] = alfa
    #             method[i] = "CenDiff"
            end
        else
                shoreAng[i] = alfa
    #             method[i] = "CenDiff"
        end
    end
                
    return shoreAng
    
end


function Snell_Law(L1,L2,alpha1)
    #     ###########################################################################    
    #     # Wave refraction using snell law.
    #     #    
    #     # INPUT:
    #     # L1:     initial wave length.
    #     # L1:     final wave length.
    #     # alpha1: initial wave dir. Cartesian notation.
    #     #
    #     # OUTPUT:
    #     # alpha1: final wave dir. Cartesian notation.
    #     ###########################################################################    
    alpha=asin.(L2.*sin.(alpha1.*pi./180.0)./L1).*180.0./pi
    
    return alpha
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

function ShoreFor(OmegaEQ,tp,hb,depthb,D50,Omega,dt,phi = 0, c = 0, D = 0, Dr = 0, Sini = 0)
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
        F = P.^.5 .* (OmegaEQ .- Omega)./std(OmegaEQ)
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

function MILLER_DEAN_CSonly(hb,depthb,sl,Y0,dt,D50,Hberm, kero, kacr, Yi)
        
    yeq = zeros(size(hb))
    Y = zeros(size(hb))

    for i in eachindex(sl)
        wl = 0.106.*hb[i].+sl[i] # # total water level
        Wast = wast(depthb[i],D50)
        yeq[i] = Y0 .- Wast.*(wl)./(Hberm.+depthb[i])
        if i ==1
            r = yeq[i].-Y[i] > 0
            k = kacr * r + kero * !r
            # k = k .* hb[i]^2
            Y[i] = Yi
        else
            r = yeq[i].-Y[i-1] > 0
            k = kacr * r + kero * !r
            # k = k .* hb[i]^2
            Y[i] = dt * k * (yeq[i] - Y[i-1]) + Y[i-1]
        end
    end        
    return Y, yeq


end

##################################################################################################
##################################################################################################
##################################################################################################


function asimiladorEnKF_MD(y,yeq,kacr,kero,kcerc,AY0,DALCS,theta,dt,DOC,Ber,Yobs,dQdx)

    #     # defining if the transect is eroding | accreting due to cross-shore forcings.
    rcross = theta.*(yeq[:,2].-y[:,2])+(1-theta).*(yeq[:,1].-y[:,1])
    kcr = zeros(size(rcross))
    kcr[rcross .>= 0] .= kacr[vec(rcross .>= 0)]
    kcr[rcross .< 0] .= kero[vec(rcross .< 0)]
    posero = rcross .< 0
    posacr = rcross .>= 0
    # kcr=kcr.*(0.5.*(hb[2:end]+hb[1:end-1]).^3)
    
    Dc=0.5.*(DOC[2:end] .+ DOC[1:end-1])

    # Yct = -kcr.*(theta.*(yeq[:,ti+1] .- y[:,ti+1]).+(1-theta).*(yeq[:,ti].- y[:,ti]))
    # Ylt = (y[:,ti+1] .- y[:,ti]).- Yct
    # dQdx = 0.5.*(dQdx[2:end]dQdx[1:end-1])
    # vlt=zeros(size(kcerc))
    Kcerc=0.5.*(kcerc[2:end] .+ kcerc[1:end-1])
    
    
    # return kalman_longitudinal_transversal(Ylt,kcerc,vlt,dQdx,dt,DALCS,Dc,Ber,sigmaK,Yct,YCTi,Yeq,kacr,kero,posero,AY0)
    return kalmanEnFK_MD(y[:,2],Kcerc,dQdx,dt,DALCS,Dc,Ber,yeq[:,2],kacr,kero,posero,posacr,AY0,Yobs)


end

function kalmanEnFK_MD(y,kcerc,dQdx,dt,DALCS,Dc,Ber,Yeq,kacr,kero,posero,posacr,AY0,Yobs)
    # function kalman_longitudinal_transversal(Ylt,kcerc,vlt,dQdx,dt,DALCS,Dc,Ber,sigmaK,Yct,YCTi,Yeq,kacr,kero,posero,AY0)
    # vectores estado
    estado_ant=[y kcerc kero.*posero kacr.*posacr AY0]' #[Ylt kcerc vlt]'
    # calculo de jacobianos
    Jacobito=calcula_jacobiano_MD(y,Dc,Ber,dQdx,DALCS["sigmaK"],kcerc,kero,kacr,DALCS["sigmaC"],Yeq,dt,posero)

    estado_post,DALCS,salto =EnFK_MD(estado_ant,Jacobito,DALCS,posero,y,Yobs)
       
    return  estado_post[1,:],estado_post[2,:],salto,DALCS,estado_post[3, posero],estado_post[4,posacr],estado_post[5,:],posero, posacr
    #Ylt,kcerc,saltoYlt,DALCS,Yct,kacr,kero,saltoYct,AY0 
    #Ylt,kcerc,vlt,saltoYlt,DALCS,Yct,kacr,kero,saltoYct,AY0
end

function calcula_jacobiano_MD(Y,dc,ber,dQdx,sigmaK,kcerc,kero,kacr,sigmaC,yeq,dt,posero)
    Jacobito=repeat(float(I(5)),1,1,length(yeq))
    # mantenemos kcerc en el jacobiano; para asegurar que no hay cambio de
    # signo
    # dF/dK1=-k0*exp(sigmak*k1)/Dc*dQdx*dt=-kcerc/Dc*dQdx*dt

    for i in eachindex(yeq)
        
        if posero[i]==1
            Jacobito[1,2,i] = -dt/(dc[i])*sigmaK[i]*kcerc[i]*dQdx[i]
            Jacobito[1,3,i] = dt*(yeq[i]-Y[i])*kero[i]*sigmaC[i]
            # Jacobito[1,3,i] = dt*(yeq[i]-Y[i])
            Jacobito[4,4,i] = 0.
            Jacobito[1,5,i] = dt*kero[i]
        else
            Jacobito[1,2,i] = -dt/(dc[i])*sigmaK[i]*kcerc[i]*dQdx[i]
            Jacobito[1,4,i] = dt*(yeq[i]-Y[i])*kacr[i]*sigmaC[i]
            # Jacobito[1,4,i] = dt*(yeq[i]-Y[i])
            Jacobito[3,3,i] = 0.
            Jacobito[1,5,i] = dt*kacr[i]
        end
    end
    return Jacobito
end

function EnFK_MD(estado_ant,Jacobito,DALCS,posero,y,Yobs)
    # Paso 1: actualizar matriz P
    H=[1 0 0 0 0]
    
    Yobs_Err=0.25
    # Hero=[1 0 1 0 0 0 0]
    # Hacr=[1 0 0 0 1 0 0]
    # nel=size(estado_ant,1)
    salto=zeros(1,size(Jacobito,3))
    estado_post=zeros(size(estado_ant))
    # estado_post[2,:] = kcerc
    # Yctn=zeros(size(Yct))
    for i in eachindex(Yobs)
        if posero[i]==1
            P0 = DALCS["P0"][i] #ensambla_lc[DALCS["P0"][i,tti],DALCS["Pero0_c"][i,tti]]
            Q = DALCS["Q"][i] #.*rand(5,5) #ensambla_lc[DALCS["Q"][i,tti],DALCS["Qero_c"][i,tti]]
            # P=cov(estado_ant;dims=2)
            P = Jacobito[:,:,i]*P0*Jacobito[:,:,i]'+ Q
             # actualizamos matrices
            DALCS["P0"][i]= P 
            HM = [1 1 1 0 1]
            # HM = [1 1 1 1 0 1 1]
        else
            P0 = DALCS["P0"][i] #ensambla_lc[DALCS["P0"][i,tti],DALCS["Pero0_c"][i,tti]]
            Q = DALCS["Q"][i] #.*rand(5,5) #ensambla_lc[DALCS["Q"][i,tti],DALCS["Qero_c"][i,tti]]
            # P=cov(estado_ant;dims=2)
            P = Jacobito[:,:,i]*P0*Jacobito[:,:,i]'+ Q
             # actualizamos matrices
            DALCS["P0"][i]= P
            HM = [1 1 0 1 1]
        end  
        # verificamos si tenemos que asimilar  
        if ! isnan(Yobs[i])
            # contador = DALCS["pos_c"][i]
            # calculamos ganancia Kalman
            # aux=H*P*H'
            R = 1 ./(length(Yobs)-1) * Yobs_Err
            K = P*H'/(H*P*H' .+ R)


            # modificamos el término en K para guardar signo
            estado_ant[2,i]=log(estado_ant[2,i]/DALCS["kcerc0"][i])/DALCS["sigmaK"][i]
            if posero[i]==1
                estado_ant[3,i]=log(estado_ant[3,i]/DALCS["kero0"][i])/DALCS["sigmaC"][i]
            else
                estado_ant[4,i]=log(estado_ant[4,i]/DALCS["kacr0"][i])/DALCS["sigmaC"][i]
            end

            modificacionKalman = K*(Yobs[i].-H*estado_ant[:,i])
            # modificacion signo de K
            estado_post[:,i] = HM'.*(estado_ant[:,i]+modificacionKalman)
            estado_post[2,i]=DALCS["kcerc0"][i]*exp(DALCS["sigmaK"][i]*estado_post[2,i])
            
            if posero[i]==1
                estado_post[3,i]=DALCS["kero0"][i]*exp(DALCS["sigmaC"][i]*estado_post[3,i])
            else
                estado_post[4,i]=DALCS["kacr0"][i]*exp(DALCS["sigmaC"][i]*estado_post[4,i])
            end
                
            
            Pnew=(float(I(length(H)))-K*H)*P
            DALCS["P0"][i]=Pnew
            salto[i]=estado_post[1,i]-y[i]
            
        else
            estado_post[:,i] = HM'.*estado_ant[:,i]
        end


    end

    return estado_post,DALCS,salto


end

##################################################################################################
##################################################################################################
##################################################################################################

function asimiladorEnKF_SF(y,kcerc,cacr,DALCS,dt,DOC,Yobs,dQdx,F,cero)

    Dc=0.5.*(DOC[2:end] .+ DOC[1:end-1])

    Kcerc=0.5.*(kcerc[2:end] .+ kcerc[1:end-1])
    
    return kalmanEnFK_SF(y,Kcerc,cacr,dQdx,dt,DALCS,Dc,Yobs,F,cero)

end

function kalmanEnFK_SF(y,kcerc,cacr,dQdx,dt,DALCS,Dc,Yobs,F,cero)
    # function kalman_longitudinal_transversal(Ylt,kcerc,vlt,dQdx,dt,DALCS,Dc,Ber,sigmaK,Yct,YCTi,Yeq,kacr,kero,posero,AY0)
    # vectores estado
    posero = F .< 0
    posacr = F .>= 0
    estado_ant = [y kcerc cero cacr]' #[Ylt kcerc vlt]'
    # calculo de jacobianos
    Jacobito=calcula_jacobiano_SF(y,Dc,dQdx,kcerc,cacr,DALCS["sigmaK"],DALCS["sigmaC"],dt,F,cero)

    estado_post,DALCS,salto =EnFK_SF(estado_ant,Jacobito,DALCS,y,Yobs,F)
       
    return  estado_post[1,:],estado_post[2,:],estado_post[3,posero],estado_post[4,posacr],salto,DALCS, posero, posacr
    #Ylt,kcerc,saltoYlt,DALCS,Yct,kacr,kero,saltoYct,AY0 
    #Ylt,kcerc,vlt,saltoYlt,DALCS,Yct,kacr,kero,saltoYct,AY0
end

function calcula_jacobiano_SF(Y,dc,dQdx,kcerc,cacr,sigmaK,sigmaC,dt,F,cero)
    Jacobito=repeat(float(I(4)),1,1,length(Y))
    # mantenemos kcerc en el jacobiano; para asegurar que no hay cambio de
    # signo
    # dF/dK1=-k0*exp(sigmak*k1)/Dc*dQdx*dt=-kcerc/Dc*dQdx*dt

    for i in eachindex(Y)
        Jacobito[1,2,i] = -dt/(dc[i])*sigmaK[i]*dQdx[i]*kcerc[i]
        if F[i]<0
            Jacobito[1,3,i] = dt*sigmaC[i]*cero[i]*(F[i])
            Jacobito[1,4,i] = 0
        else
            Jacobito[1,3,i] = 0
            Jacobito[1,4,i] = dt*sigmaC[i]*cacr[i]*(F[i])
        end
    end
    return Jacobito
end

function EnFK_SF(estado_ant,Jacobito,DALCS,y,Yobs,F)
    # Paso 1: actualizar matriz P
    H=[1 0 0 0]
    Yobs_Err=0.25 #.*rand(length(Yobs),1)
    salto=zeros(1,size(Jacobito,3))
    estado_post=zeros(size(estado_ant))

    for i in eachindex(Yobs)

        P0 = DALCS["P0"][i] #ensambla_lc[DALCS["P0"][i,tti],DALCS["Pero0_c"][i,tti]]
        Q = DALCS["Q"][i] #.*rand(3,3) #ensambla_lc[DALCS["Q"][i,tti],DALCS["Qero_c"][i,tti]]
        P = Jacobito[:,:,i]*P0*Jacobito[:,:,i]'+Q
            # actualizamos matrices
        DALCS["P0"][i]= P
        if F[i]<0
            HM = [1 1 1 0]
        else
            HM = [1 1 0 1]
        end
        
        if ! isnan(Yobs[i])

            R = 1 ./(length(Yobs)-1) * Yobs_Err #' * Yobs_Err
            K = P*H'/(H*P*H' .+ R)

            estado_ant[2,i]=log(estado_ant[2,i]/DALCS["kcerc0"][i])/DALCS["sigmaK"][i]
            if F[i]<0
                estado_ant[3,i]=log(estado_ant[3,i]/DALCS["c0"][i])/DALCS["sigmaC"][i]
            else
                estado_ant[4,i]=log(estado_ant[4,i]/DALCS["c0"][i])/DALCS["sigmaC"][i]
            end

            modificacionKalman = K*(Yobs[i].-H*estado_ant[:,i])
            # modificacion signo de K
            estado_post[:,i] = HM'.*(estado_ant[:,i]+modificacionKalman)
            estado_post[2,i]=DALCS["kcerc0"][i]*exp(DALCS["sigmaK"][i]*estado_post[2,i])
            if F[i]<0
                estado_post[3,i]=DALCS["c0"][i]*exp(DALCS["sigmaC"][i]*estado_post[3,i])
            else
                estado_post[4,i]=DALCS["c0"][i]*exp(DALCS["sigmaC"][i]*estado_post[4,i])
            end
            
            Pnew=(float(I(length(H)))-K*H)*P
            DALCS["P0"][i]=Pnew
            salto[i]=estado_post[1,i]-y[i]
            
        else
            estado_post[:,i] = HM'.*estado_ant[:,i]
        end


    end

    return estado_post,DALCS,salto


end

##################################################################################################
##################################################################################################
##################################################################################################


function asimiladorEnKF_Y09(y,hs,kcerc,HSB,DY,DT,DALCS,dt,DOC,Yobs,dQdx)

    Dc=0.5.*(DOC[2:end] .+ DOC[1:end-1])

    # Yct = -kcr.*(theta.*(yeq[:,ti+1] .- y[:,ti+1]).+(1-theta).*(yeq[:,ti].- y[:,ti]))
    # Ylt = (y[:,ti+1] .- y[:,ti]).- Yct
    # dQdx = 0.5.*(dQdx[2:end]dQdx[1:end-1])
    # vlt=zeros(size(kcerc))
    Kcerc = 0.5.*(kcerc[2:end] .+ kcerc[1:end-1])
    hss = 0.5.*(hs[2:end] .+ hs[1:end-1])
    # hss[1] = hss[2] 
    # hss[end] = hss[end-1]
    
    # return kalman_longitudinal_transversal(Ylt,kcerc,vlt,dQdx,dt,DALCS,Dc,Ber,sigmaK,Yct,YCTi,Yeq,kacr,kero,posero,AY0)
    return kalmanEnFK_Y09(y,hss,Kcerc,HSB,DY,DT,dQdx,dt,DALCS,Dc,Yobs)


end

function kalmanEnFK_Y09(y,hs,kcerc,HSB,DY,DT,dQdx,dt,DALCS,Dc,Yobs)
    # function kalman_longitudinal_transversal(Ylt,kcerc,vlt,dQdx,dt,DALCS,Dc,Ber,sigmaK,Yct,YCTi,Yeq,kacr,kero,posero,AY0)
    # vectores estado
    estado_ant=[y kcerc HSB DT DY]' #[Ylt kcerc vlt]'
    # calculo de jacobianos
    Jacobito=calcula_jacobiano_Y09(y,hs,Dc,dQdx,kcerc,HSB,DY,DT,DALCS["sigmaK"],DALCS["sigmaHSB"],
    DALCS["sigmaDT"],DALCS["sigmaDY"],dt)

    estado_post,DALCS,salto =EnFK_Y09(estado_ant,Jacobito,DALCS,y,Yobs)
       
    return  estado_post[1,:],estado_post[2,:],estado_post[3, :],estado_post[4, :],estado_post[5,:],salto,DALCS
    #Ylt,kcerc,saltoYlt,DALCS,Yct,kacr,kero,saltoYct,AY0 
    #Ylt,kcerc,vlt,saltoYlt,DALCS,Yct,kacr,kero,saltoYct,AY0
end

function calcula_jacobiano_Y09(Y,hs,dc,dQdx,kcerc,HSB,DY,DT,sigmaK,sigmaHSB,sigmaDT,sigmaDY,dt)
    Jacobito=repeat(float(I(5)),1,1,length(Y))
    # mantenemos kcerc en el jacobiano; para asegurar que no hay cambio de
    # signo
    # dF/dK1=-k0*exp(sigmak*k1)/Dc*dQdx*dt=-kcerc/Dc*dQdx*dt

    for i in eachindex(Y)
        Jacobito[1,2,i] = -dt/(dc[i])*sigmaK[i]*kcerc[i]*dQdx[i]
        Jacobito[1,3,i] = dt/(DT[i]*HSB[i])*sigmaHSB[i]*hs[i]*(3*hs[i]^2*DY[i]/HSB[i]^2+Y[i]+DY[i])
        Jacobito[1,4,i] = dt/DT[i]*sigmaDT[i]*hs[i]/HSB[i]*(DY[i]*(hs[i]^2-HSB[i]^2)/HSB[i]^2+Y[i]) 
        Jacobito[1,5,i] = dt*DY[i]*sigmaDY[i]/DT[i]*(hs[i]^3-HSB[i]^2*hs[i])/(HSB[i]^3)
    end
    return Jacobito
end

function EnFK_Y09(estado_ant,Jacobito,DALCS,y,Yobs)
    # Paso 1: actualizar matriz P
    H=[1 0 0 0 0]
    Yobs_Err=0.25 #.*rand(length(Yobs),1)
    # Hero=[1 0 1 0 0 0 0]
    # Hacr=[1 0 0 0 1 0 0]
    # nel=size(estado_ant,1)
    salto=zeros(1,size(Jacobito,3))
    estado_post=zeros(size(estado_ant))
    # estado_post[2,:] = kcerc
    # Yctn=zeros(size(Yct))
    for i in eachindex(Yobs)
        P0 = DALCS["P0"][i] #ensambla_lc[DALCS["P0"][i,tti],DALCS["Pero0_c"][i,tti]]
        Q = DALCS["Q"][i] #.*rand(5,5) #ensambla_lc[DALCS["Q"][i,tti],DALCS["Qero_c"][i,tti]]
        P = Jacobito[:,:,i]*P0*Jacobito[:,:,i]'+Q
        DALCS["P0"][i]= P 
        HM = [1 1 1 1 1]
        # verificamos si tenemos que asimilar  
        if ! isnan(Yobs[i])
            # contador = DALCS["pos_c"][i]
            # calculamos ganancia Kalman
            # aux=H*P*H'
            R = 1 ./(length(Yobs)-1) * Yobs_Err #' * Yobs_Err
            K = P*H'/(H*P*H'.+ R)
            # modificamos el término en K para guardar signo
            estado_ant[2,i]=log(estado_ant[2,i]/DALCS["kcerc0"][i])/DALCS["sigmaK"][i]
            estado_ant[3,i]=log(estado_ant[3,i]/DALCS["HSB0"][i])/DALCS["sigmaHSB"][i]
            estado_ant[4,i]=log(estado_ant[4,i]/DALCS["DT0"][i])/DALCS["sigmaDT"][i]
            estado_ant[5,i]=log(estado_ant[5,i]/DALCS["DY0"][i])/DALCS["sigmaDY"][i]

            modificacionKalman = K*(Yobs[i].-H*estado_ant[:,i])
            # modificacion signo de K
            estado_post[:,i] = HM'.*(estado_ant[:,i]+modificacionKalman)
            estado_post[2,i]=DALCS["kcerc0"][i]*exp(DALCS["sigmaK"][i]*estado_post[2,i])
            estado_post[3,i]=DALCS["HSB0"][i]*exp(DALCS["sigmaHSB"][i]*estado_post[3,i])
            estado_post[4,i]=DALCS["DT0"][i]*exp(DALCS["sigmaDT"][i]*estado_post[4,i])
            estado_post[5,i]=DALCS["DY0"][i]*exp(DALCS["sigmaDY"][i]*estado_post[5,i])
            
            # if posero[i]==1
            #     estado_post[3,i]=DALCS["kero0"][i]*exp(DALCS["sigmaC"][i]*estado_post[3,i])
            # else
            #     estado_post[4,i]=DALCS["kacr0"][i]*exp(DALCS["sigmaC"][i]*estado_post[4,i])
            # end
                
            
            Pnew=(float(I(length(H)))-K*H)*P
            DALCS["P0"][i]=Pnew
            salto[i]=estado_post[1,i]-y[i]
            
        else
            estado_post[:,i] = HM'.*estado_ant[:,i]
        end


    end

    return estado_post,DALCS,salto


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