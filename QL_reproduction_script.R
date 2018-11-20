rm(list = ls())
                                        
sim_name <- "QL_reproduction"

## get and create directories, load functions file
working_directory <- getwd()
bash_date <- format(Sys.time(),format="%Y%m%d")
function_file <- sprintf("%s/%s_%s_functions.R", working_directory, bash_date, sim_name)
dir.create(file.path(working_directory, "/user_data/"), showWarnings = FALSE)
output_directory <- sprintf("%s/user_data/", getwd())
source(function_file)

## LOAD PACKAGES ##
require("doMC")
require("foreach")

## PARAMETERS ##
## number of alleles in population
number_alleles <- 15
## number of drones each queen mates with
number_drone_matings <- 25
## cost to colony of producing diploid (dead) males
cost_homozoygosity <- 0.05
## number of generations (each generation = 6 months)
number_generations <- 6
## proportion of swarms that become queen-less
probability_QL_colony <- 0.3
## probability that queen AND colony survives each generation
prob_queen_survives <- 1 - probability_QL_colony
## mean swarms per colony
average_swarms <- 4
## proportion of production of drones by QL colonies (relative to QR)
QL_drone_production <- 1
## initial distribution alleles (assume alleles are at balancing selection equilibrium - i.e. even frequency)
initial_distribution_alleles <- rep(1 / number_alleles,number_alleles) 
## number of monte carlo replicates
num_trials <- 100
## number of cores
num_cores <- 2
## set up parallel backend
registerDoMC(cores = num_cores)

results <- foreach (loop = 1:num_trials, .combine = rbind) %dopar% { 
    
    repeat { ## if the population goes extinct, I want to repeat the simulation (to standarize the number of "successful" invasions  

        temp_list <- list()

        ## initialise objects to store information about the simulation
        queen_allele_frequencies <- matrix(nrow = 0, ncol = number_alleles)
        worker_laid_drone_allele_distribution <- matrix(nrow = 0, ncol = number_alleles)
        queen_laid_drone_allele_distribution <- matrix(nrow = 0, ncol = number_alleles)
        simulation_extinction_status <- NULL ## 0 for extinct; 1 for no extinction
        generation_all_QR_colonies_lost <- NULL
        population_size <- NULL

        ## INITIALISE INVADING COLONY ##
        
        ## produce population matrix 
        population <- matrix(numeric(2 + number_alleles + 2),nrow = 1, ncol = 2 + number_alleles + 2)

        ## choose queen alleles from the source population
        population <- initialiseQueenAllelesFromSourcePop(population, number_alleles, initialise_distribution_alleles)

        ## choose which drones the invading queen mates with (also from source population)
        population <- chooseDroneAlleles(population, colony_ID = 1, number_alleles, number_drone_matings, initial_distribution_alleles)
        
        ## calculate colony fitness of invading colony
        population <- calculateColonyFitness(population, colony_ID = 1, number_alleles, cost_homozygosity)
        
        ## first colony is queenright
        population <- setColonyQueenStatus(population, number_alleles, colony_ID = 1, queen_status = 1)

        ## record simulation data (queen alleles, worker-laid drone allels, queen-laid drone alleles, population size)
        queen_allele_frequencies <- rbind( queen_allele_frequencies,
                                          tabulate(bin = population[population[ , number_alleles + 4] == 1, 1:2], nbins = number_alleles) /
                                          ( NROW( population[population[, number_alleles + 4] == 1] ) * 2 ) )
        queen_laid_allele_distribution <- rbind( queen_laid_allele_distribution, unlist(list_output[2]) ) ## queen-laid drones
        worker_laid_allele_distribution <- rbind( worker_laid_allele_distribution, unlist(list_output[3]) ) ## worker-laid drones
        population_size[length(population_size) + 1] <- NROW(population) ## number of colonies left

        ## INVASION ##
        for (generationloop in 1:number_generations){
            
            ## logistic growth rate (depends on population size, we include both QL and QR colonies)
            average_growth_rate <- average_swarms 
                                        # growth rate weighted by colony fitness (only QR)
            weighted_growth_rate <- average_growth_rate * population[,number_alleles + 3] 
            
                                        # how many daughter swarms does each colony leave?
            swarming_vector <- rpois(n = N, lambda = weighted_growth_rate)
                                        #generate new population
            population <- generate_new_population(swarming_vector,population,N,number_alleles,number_drone_matings,
                                                  prob_queen_survives,cost_homozoygosity,probability_QL_colony)
            
                                        #is the population extinct?
            if (population == 'extinct'){ 
## i want to record data from those simulations with extinct poulations now                
                counter <- counter + 1

                break #breaks out of for 1:number_generations loop
                
            } else { #not extinct, proceed normally
                
                                        #no longer need the old population size 
                N <- nrow(population)
                
                                        #record information about the simulation
                queen_allele_frequencies[counter,] <- sort(tabulate(bin = population[,1:2],nbins = number_alleles) /
                                                           (N * 2), decreasing = TRUE) #counts, ordered from highest freq to lowest
                
                                        #drone alleles are stored as proportions not counts so use colMeans
                                        #if N == 1, colMeans throws an error (array becomes vector and it requires 2 dim object) so use drop = FALSE
                drone_allele_frequencies[counter,] <- sort(colMeans(population[,3:(number_alleles + 2), drop = FALSE]),decreasing = TRUE) 
                
                counter <- counter + 1
            } 
        }
        
                                        #this code is executed either because the break was triggered (population is extinct) or through normal 
                                        #termination of the for loop. I don't want to record the former but want to record the latter.
        
        if (population != 'extinct'){ #normal termination
                                        #break out of repeat loop
            
            break # transfers control to foreach loop
            
        } #if population == 'extinct', record nothing and go back to start of repeat loop
        
                                        #this code is only executed if the population == 'extinct'
        extinction_counter <- extinction_counter + 1
        ## record some information here
        
    } #this marks the end of the repeat loop
    
                                        #record simulation 
    
    temp_list[[1]] <- queen_allele_frequencies
    temp_list[[2]] <- drone_allele_frequencies
    
    temp_list[[4]] <- extinction_counter
    
    results <- temp_list
    
} #this marks the end of the foreach loop

filename <- sprintf("%s/%s_%s_na%d_ch%d_as%d_pq%d_dm%d.RData",
                    output_directory,bash_date,sim_name,number_alleles,cost_homozoygosity*100,average_swarms,probability_QL_colony*100,number_drone_matings)
save(results, file = filename)


