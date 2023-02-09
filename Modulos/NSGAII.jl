using Distributions
using Optim
using Plots

function nsga2(model, obj_fun, num_params, num_generations, num_individuals)
    # Initialize population of candidate solutions
    population = rand(Normal(0, 1), num_individuals, num_params)
    
    # Evaluate initial population
    objectives = obj_fun.(population)
    
    # Initialize iteration counter
    t = 1
    
    while t <= num_generations
        # Select parents for breeding
        parents = tournament_selection(population, objectives)
        
        # Create offspring using genetic operators
        offspring = genetic_operators(parents)
        
        # Evaluate offspring
        offspring_obj = obj_fun.(offspring)
        
        # Update population with offspring and original solutions
        population, objectives = non_dominated_sort(
            vcat(population, offspring),
            vcat(objectives, offspring_obj)
        )
        
        # Truncate population to desired size
        population = population[1:num_individuals, :]
        objectives = objectives[1:num_individuals, :]
        
        # Update iteration counter
        t += 1
    end
    
    return population, objectives
end

function tournament_selection(population, objectives)
    num_individuals, num_params = size(population)
    num_selected = convert(Int, num_individuals / 2)
    selected = zeros(num_selected, num_params)
    
    for i in 1:num_selected
        competitors = rand(1:num_individuals, 2, 1)
        if dominates(objectives[competitors[1], :], objectives[competitors[2], :])
            selected[i, :] = population[competitors[1], :]
        else
            selected[i, :] = population[competitors[2], :]
        end
    end
    
    return selected
end

function genetic_operators(parents)
    num_parents, num_params = size(parents)
    num_offspring = num_parents
    offspring = zeros(num_offspring, num_params)
    
    for i in 1:2:num_parents
        # Select parents to breed
        p1 = parents[i, :]
        p2 = parents[i + 1, :]
        
        # Apply crossover operator
        o1, o2 = crossover(p1, p2)
        
        # Apply mutation operator
        o1 = mutation(o1)
        o2 = mutation(o2)
        
        offspring[i, :] = o1
        offspring[i + 1, :] = o2
    end
    
    return offspring
end

function crossover(p1, p2)
    alpha = rand()
    o1 = alpha * p1 + (1 - alpha) * p2
    o2 = (1 - alpha) * p1 + alpha * p2
    return o1, o2
end

function dominates(obj1, obj2)
    # Determine if obj1 dominates obj2
    # A solution dominates another if it is better in all objectives
    # and not worse in at least one objective.
    
    dominates = true
    worse = false
    
    for i in 1:length(obj1)
        if obj1[i] > obj2[i]
            dominates = false
            break
        elseif obj1[i] < obj2[i]
            worse = true
        end
    end
    
    return dominates && worse
end

function non_dominated_sort(population, objectives)
    # Sort population into non-dominated fronts
    fronts = sort_fronts(population, objectives)
    
    # Select individuals from each front
    selected = select_individuals(fronts)
    
    # Update population and objectives arrays
    population = population[selected, :]
    objectives = objectives[selected, :]
    
    return population, objectives
end

function sort_fronts(population, objectives)
    # Sort population into non-dominated fronts
    num_individuals, num_obj = size(objectives)
    
    fronts = []
    front = []
    
    for p in 1:num_individuals
        dominated = false
        
        for q in front
            if dominates(objectives[q, :], objectives[p, :])
                dominated = true
                break
            end
        end
        
        if !dominated
            for q in 1:num_individuals
                if p == q
                    continue
                end
                
                if dominates(objectives[p, :], objectives[q, :])
                    deleteat!(front, findfirst(front .== q))
                end
            end
            
            push!(front, p)
        end
        
        if length(front) + length(fronts) == num_individuals
            break
        end
        
        if length(front) == 0
            push!(fronts, front)
            front = []
        end
    end
    
    push!(fronts, front)
    
    return fronts
end

function select_individuals(fronts)
    # Select individuals from each front
    selected = []
    
    for front in fronts
        n = length(front)
        if n > 0
            if n <= length(selected)
                selected[1:n] = front
            else
                selected = vcat(selected, front)
            end
        end
    end
    
    return selected
end

function mutation(o)
    # Apply mutation operator to offspring
    # Add Gaussian noise with mean 0 and standard deviation 1 to each gene
    mutated = o .+ randn(size(o))
    return mutated
end
