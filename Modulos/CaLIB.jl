module CAL

using Polynomials
using LinearAlgebra
using SciPy
using Statistics
using ImageFiltering
using NumericalIntegration
using Shuffle



# Here are the steps to calibrate a model using genetic algorithms:

# 1. Define the model parameters that need to be optimized.
# 2. Define a fitness function that measures the quality of the model predictions given a set of parameters.
# 3. Initialize a population of parameter sets, either randomly or using some prior knowledge.
# 4. Evaluate the fitness of each member of the population using the fitness function.
# 5. Select the best individuals from the population based on their fitness scores.
# 6. Breed the selected individuals to produce a new generation of parameter sets, using genetic operations such as crossover and mutation.
# 7. Repeat steps 4 to 6 for a specified number of generations or until the fitness function converges to a satisfactory solution.
# 8. Return the best set of parameters found.
# This process is a heuristic optimization method and can be time-consuming,
# but it can be very effective for complex, non-linear models with many parameters.
# It's also a good option when traditional optimization methods, such as gradient descent,
# struggle to find a good solution.


function genetic_algorithm(fitness_function, parameter_bounds, population_size, 
    max_generations, mutation_rate)
    population = [rand(parameter_bounds[i][1]:(parameter_bounds[i][2] - parameter_bounds[i][1]), population_size)
                  for i in 1:length(parameter_bounds)]
    for generation in 1:max_generations
        fitness_values = fitness_function.(population)
        sorted_population = sort(zip(fitness_values, population), by=x -> x[1], rev=true)
        elites = [x[2] for x in sorted_population[1:round(Int, population_size / 2)]]
        new_population = []
        while length(new_population) < population_size
            parent1, parent2 = elites[rand(1:length(elites))], elites[rand(1:length(elites))]
            child = [0.5 * (parent1[i] + parent2[i]) + mutation_rate * randn()
                     for i in 1:length(parameter_bounds)]
            if all(x -> x >= parameter_bounds[i][1] && x <= parameter_bounds[i][2], child)
                push!(new_population, child)
            end
        end
        population = new_population
    end
    return sorted_population[1][2]
end

# The S
# This function takes a 'fitness_function', which is the function to be optimized, 
# a 'parameter_bounds' array that specifies the bounds for each parameter,
# 'population_size', which is the number of parameter sets in each generation, 
# 'max_generations' which is the maximum number of generations to run the 
# algorithm for, and 'mutation_rate', which determines the magnitude of the random
# perturbations applied during mutation. The function returns the best set of
# parameters found by the genetic algorithm.


function sce_ua(f, x0, ngen, npop, npar)
    # Initialize the population
    pop = Array{Float64}(undef, npop, npar)
    for i in 1:npop
        pop[i,:] = x0 + randn(npar)
    end
    
    # Main loop
    for gen in 1:ngen
        # Shuffle the population
        shuffle!(pop)
        
        # Complex evolution
        for i in 1:2:npop
            x1 = pop[i,:]
            x2 = pop[i+1,:]
            y1 = f(x1)
            y2 = f(x2)
            
            if y1 < y2
                pop[i+1,:] = x1 + 0.5 * (x2 - x1) + 0.5 * randn(npar)
            else
                pop[i,:] = x2 + 0.5 * (x1 - x2) + 0.5 * randn(npar)
            end
        end
    end
    
    # Return the best solution
    return pop[argmin(f.(pop))[1],:]
end

# SCE-UA algorithm in Julia:

# f: This is the function that you want to optimize. It should take a single input argument x (of length npar) and return a scalar value.

# x0: This is the initial guess for the solution to the optimization problem. It should be a 1-dimensional array of length npar.

# ngen: This is the number of generations (i.e., iterations) of the SCE-UA algorithm.

# npop: This is the number of individuals (candidate solutions) in the population. It should be an even number, as the SCE-UA algorithm requires that the population size be divisible by 2.

# npar: This is the number of parameters in the optimization problem. It should be the length of x0 and the number of elements in each row of the pop array.


function sce_ua2(f, x0, ngen, npop, npar, mag)
    # Initialize the population
    pop = Array{Float64}(undef, npop, npar)
    for i in 1:npop
        pop[i,:] = x0 .+  mag .* (2. * rand(npar) .- 1)
        # for j in 1:npar
        #     pop[i,j] = min(upper_bounds[j], max(lower_bounds[j], pop[i,j]))
        # end
    end


    # Main loop
    fvals = Array{Float64}(undef, npop)
    for gen in 1:ngen
        # Shuffle the population
        for i in 1:npar
            pop[:,i] = shuffle(pop[:,i])
        end
        
        # Evaluate the function for all individuals in the population
        for i in 1:npop
            fvals[i] = f(pop[i,:])
        end
        
        # Complex evolution
        for i in 1:2:npop
            x1 = pop[i,:]
            x2 = pop[i+1,:]
            y1 = fvals[i]
            y2 = fvals[i+1]
            
            if y1 < y2
                # Use x1 as the better individual
                pop[i+1,:] = x1 .+ 0.5 .* (x2 .- x1) .+ 0.5 .* mag .* (2. * rand(npar) .- 1)
            else
                # Use x2 as the better individual
                pop[i,:] = x2 .+ 0.5 .* (x1 .- x2) .+ 0.5 .* mag .* (2. * rand(npar) .- 1)
            end
        end
        println("Progress = " * string(gen/ngen*100) * " %")
        println("Generation: $gen / $ngen")
        println()
    end
    
    # Return the best solution
    return pop[argmin(fvals)[1],:], minimum(fvals)
end

####################################################
################## Kalman Filters ##################
####################################################


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
end # module