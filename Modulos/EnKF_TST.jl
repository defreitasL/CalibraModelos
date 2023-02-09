
module EnKF
##################################################################################################
##################################################################################################
##################################################################################################
using Statistics

function asimiladorEnKF_MD(y,yeq,kacr,kero,kcerc,AY0,DALCS,theta,Yobs,err_mag)

    #     # defining if the transect is eroding | accreting due to cross-shore forcings.
    rcross = theta.*(yeq[:,2].-y[:,2])+(1-theta).*(yeq[:,1].-y[:,1])
    kcr = zeros(size(rcross))
    kcr[rcross .>= 0] .= kacr[vec(rcross .>= 0)]
    kcr[rcross .< 0] .= kero[vec(rcross .< 0)]
    posero = rcross .< 0
    posacr = rcross .>= 0
    # kcr=kcr.*(0.5.*(hb[2:end]+hb[1:end-1]).^3)
    
    # Dc=0.5.*(DOC[2:end] .+ DOC[1:end-1])

    # Yct = -kcr.*(theta.*(yeq[:,ti+1] .- y[:,ti+1]).+(1-theta).*(yeq[:,ti].- y[:,ti]))
    # Ylt = (y[:,ti+1] .- y[:,ti]).- Yct
    # dQdx = 0.5.*(dQdx[2:end]dQdx[1:end-1])
    # vlt=zeros(size(kcerc))
    Kcerc=0.5.*(kcerc[2:end] .+ kcerc[1:end-1])
    
    
    # return kalman_longitudinal_transversal(Ylt,kcerc,vlt,dQdx,dt,DALCS,Dc,Ber,sigmaK,Yct,YCTi,Yeq,kacr,kero,posero,AY0)
    return kalmanEnFK_MD(y[:,2],Kcerc,DALCS,kacr,kero,posero,posacr,AY0,Yobs,err_mag)


end

function kalmanEnFK_MD(y,kcerc,DALCS,kacr,kero,posero,posacr,AY0,Yobs,err_mag)
    # function kalman_longitudinal_transversal(Ylt,kcerc,vlt,dQdx,dt,DALCS,Dc,Ber,sigmaK,Yct,YCTi,Yeq,kacr,kero,posero,AY0)
    # vectores estado
    kcerc = 1.0 ./ vec(DALCS["sigmaK"]') .* log.(kcerc./vec(DALCS["kcerc0"]'))
    kacr = 1.0 ./ vec(DALCS["sigmaC"]') .* log.(kacr./vec(DALCS["kacr0"]'))
    kero = 1.0 ./ vec(DALCS["sigmaC"]') .* log.(kero./vec(DALCS["kero0"]'))

    estado_ant=[y kcerc kero.*posero kacr.*posacr AY0]' #[Ylt kcerc vlt]'
    # calculo de jacobianos
    # Jacobito=calcula_jacobiano_MD(y,Dc,Ber,dQdx,DALCS["sigmaK"],kcerc,kero,kacr,DALCS["sigmaC"],Yeq,dt,posero)

    estado_post, salto =EnFK_MD(estado_ant,DALCS,Yobs,err_mag,posero,posacr)
       
    return  estado_post[1,:],estado_post[2,:],estado_post[3, posero],estado_post[4,posacr],estado_post[5,:],posero, posacr,salto,DALCS
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
            Jacobito[1,2,i] = -dt/(dc[i])*sigmaK[i]*kcerc[i]*dQdx[i]-kero[i]
            Jacobito[1,3,i] = dt*(yeq[i]-Y[i])#*kero[i]*sigmaC[i]
            # Jacobito[1,3,i] = dt*(yeq[i]-Y[i])
            Jacobito[4,4,i] = 0.
            Jacobito[1,5,i] = dt*kero[i]
        else
            Jacobito[1,2,i] = -dt/(dc[i]+ber[i])*sigmaK[i]*kcerc[i]*dQdx[i]-kacr[i]
            Jacobito[1,4,i] = dt*(yeq[i]-Y[i])#*kacr[i]*sigmaC[i]
            # Jacobito[1,4,i] = dt*(yeq[i]-Y[i])
            Jacobito[3,3,i] = 0.
            Jacobito[1,5,i] = dt*kacr[i]
        end
    end
    return Jacobito
end

function EnFK_MD(estado_ant,DALCS,Yobs,err_mag,posero,posacr)
    # Paso 1: actualizar matriz P
    H=[1 0 0 0 0]
    # if posero[i]==1
    #     HM = [1 1 1 0 1]
    # else
    #     HM = [1 1 0 1 1]
    # end
    salto=zeros(size(Yobs))
    estado_post=zeros(size(estado_ant))
    P=zeros(size(H))
    Yobs_err = err_mag.*rand(length(Yobs))
    # estado_post[2,:] = kcerc
    # Yctn=zeros(size(Yct))
    P=cov(estado_ant;dims=2)
    R1=1/(length(Yobs)-1)*Yobs_err'*Yobs_err
    Kal=(P*H')/(H*P*H'.+R1)
    for i in eachindex(Yobs)
        Ymod = estado_ant[:,i].*H'
        sl_err = Yobs[i]+Yobs_err[i]-Ymod[1]
        if posero[i]==1
            HM = [1 1 1 0 1]
        else
            HM = [1 1 0 1 1]
        end
        estado_post[:,i] .= (estado_ant[:,i].+ Kal .* sl_err) .* HM'
    end
    estado_post[2,:] .= vec(DALCS["kcerc0"]).*exp.(vec(DALCS["sigmaK"]).*estado_post[2,:])
    estado_post[3,:] .= vec(DALCS["kero0"]).*exp.(vec(DALCS["sigmaC"]).*estado_post[3,:])
    estado_post[4,:] .= vec(DALCS["kacr0"]).*exp.(vec(DALCS["sigmaC"]).*estado_post[4,:])
    salto .= estado_post[1,:] .- estado_ant[1,:]

    return estado_post,salto


end

##################################################################################################
##################################################################################################
##################################################################################################

function asimiladorEnKF_SF(y,kcerc,c,DALCS,dt,DOC,Yobs,dQdx,F,r)

    Dc=0.5.*(DOC[2:end] .+ DOC[1:end-1])
    Dc[1]=Dc[2]
    Dc[end]=Dc[end-1]
    Kcerc=0.5.*(kcerc[2:end] .+ kcerc[1:end-1])
    kcerc[1]=kcerc[2]
    kcerc[end]=kcerc[end-1]
    return kalmanEnFK_SF(y,Kcerc,c,dQdx,dt,DALCS,Dc,Yobs,F,r)

end

function kalmanEnFK_SF(y,kcerc,c,dQdx,dt,DALCS,Dc,Yobs,F,r)
    # function kalman_longitudinal_transversal(Ylt,kcerc,vlt,dQdx,dt,DALCS,Dc,Ber,sigmaK,Yct,YCTi,Yeq,kacr,kero,posero,AY0)
    # vectores estado
    estado_ant=[y kcerc c]' #[Ylt kcerc vlt]'
    # calculo de jacobianos
    Jacobito=calcula_jacobiano_SF(y,Dc,dQdx,kcerc,c,DALCS["sigmaK"],DALCS["sigmaC"],dt,F,r)

    estado_post,DALCS,salto =EnFK_SF(estado_ant,Jacobito,DALCS,y,Yobs,F)
       
    return  estado_post[1,:],estado_post[2,:],estado_post[3,:],salto,DALCS
    #Ylt,kcerc,saltoYlt,DALCS,Yct,kacr,kero,saltoYct,AY0 
    #Ylt,kcerc,vlt,saltoYlt,DALCS,Yct,kacr,kero,saltoYct,AY0
end

function calcula_jacobiano_SF(Y,dc,dQdx,kcerc,c,sigmaK,sigmaC,dt,F,r)
    Jacobito=repeat(float(I(3)),1,1,length(Y))
    # mantenemos kcerc en el jacobiano; para asegurar que no hay cambio de
    # signo
    # dF/dK1=-k0*exp(sigmak*k1)/Dc*dQdx*dt=-kcerc/Dc*dQdx*dt

    for i in eachindex(Y)
        Jacobito[1,2,i] = -dt/(dc[i])*sigmaK[i]*dQdx[i]*kcerc[i]
        if F[i]<0
            Jacobito[1,3,i] = dt*sigmaC[i]*c[i]*(r[i]*F[i])
        else
            Jacobito[1,3,i] = dt*sigmaC[i]*c[i]*(F[i])
        end
    end
    return Jacobito
end

function EnFK_SF(estado_ant,Jacobito,DALCS,y,Yobs,F)
    # Paso 1: actualizar matriz P
    H=[1 0 0]
    salto=zeros(1,size(Jacobito,3))
    estado_post=zeros(size(estado_ant))

    for i in eachindex(Yobs)

        P0 = DALCS["P0"][i] #ensambla_lc[DALCS["P0"][i,tti],DALCS["Pero0_c"][i,tti]]
        Q = DALCS["Q"][i] #ensambla_lc[DALCS["Q"][i,tti],DALCS["Qero_c"][i,tti]]
        P = Jacobito[:,:,i]*P0*Jacobito[:,:,i]'+Q
            # actualizamos matrices
        DALCS["P0"][i]= P 
        HM = [1 1 1]
        # HM = [1 1 1 1 0 1 1]
        
        if ! isnan(Yobs[i])

            K = P*H'/(H*P*H'.+ DALCS["R"][i])

            estado_ant[2,i]=log(estado_ant[2,i]/DALCS["kcerc0"][i])/DALCS["sigmaK"][i]
            estado_ant[3,i]=log(estado_ant[3,i]/DALCS["c0"][i])/DALCS["sigmaC"][i]

            modificacionKalman = K*(Yobs[i].-H*estado_ant[:,i])
            # modificacion signo de K
            estado_post[:,i] = HM'.*(estado_ant[:,i]+modificacionKalman)
            estado_post[2,i]=DALCS["kcerc0"][i]*exp(DALCS["sigmaK"][i]*estado_post[2,i])
            estado_post[3,i]=DALCS["c0"][i]*exp(DALCS["sigmaC"][i]*estado_post[3,i])
            
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
    hss[1] = hss[2] 
    hss[end] = hss[end-1]
    
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

    # Hero=[1 0 1 0 0 0 0]
    # Hacr=[1 0 0 0 1 0 0]
    # nel=size(estado_ant,1)
    salto=zeros(1,size(Jacobito,3))
    estado_post=zeros(size(estado_ant))
    # estado_post[2,:] = kcerc
    # Yctn=zeros(size(Yct))
    for i in eachindex(Yobs)
        P0 = DALCS["P0"][i] #ensambla_lc[DALCS["P0"][i,tti],DALCS["Pero0_c"][i,tti]]
        Q = DALCS["Q"][i] #ensambla_lc[DALCS["Q"][i,tti],DALCS["Qero_c"][i,tti]]
        P = Jacobito[:,:,i]*P0*Jacobito[:,:,i]'+Q
        DALCS["P0"][i]= P 
        HM = [1 1 1 1 1]
        # verificamos si tenemos que asimilar  
        if ! isnan(Yobs[i])
            # contador = DALCS["pos_c"][i]
            # calculamos ganancia Kalman
            # aux=H*P*H'
            K = P*H'/(H*P*H'.+ DALCS["R"][i])

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


end #module