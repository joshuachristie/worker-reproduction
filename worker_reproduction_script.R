rm(list = ls())

## LOAD PACKAGES ##
require("doMC")
require("foreach")

################
## PARAMETERS ##
################

## MODEL PARAMETERS ##
## number of alleles in population
number_alleles <- 15
## number of drones each queen mates with
number_drone_matings <- 25
## cost to colony of producing diploid (dead) males
cost_homozygosity <- 0.05
#proportion of swarms that become queen-less
probability_QL_colony <- 0.3
#probability that queen AND colony survives each generation
prob_queen_survives <- 1 - probability_QL_colony
## proportion of swarms that become queen-less
probability_QL_colony <- 0.3
## probability that queen AND colony survives each generation
prob_queen_survives <- 1 - probability_QL_colony
## mean swarms per colony
average_swarms <- 2
## worker reproduction status (TRUE = QL colonies produce workers; FALSE = QL colonies do not produce workers)
worker_reproduction_status <- TRUE
## proportion of production of drones by QL colonies (relative to QR)
QL_drone_production <- 1
## distribution alleles in source population (we assume alleles in source population are under strong balancing selection)
initial_distribution_alleles <- rep(1 / number_alleles,number_alleles) ## each allele present at the same frequency
## are there multiple founders or a single founder? (if the former, a second colony invades after a delay)
multiple_founder_status <- TRUE
## total number of generations (each generation = 6 months)
number_generations <- 6
## number of generations before second invasion (if applicable)
generations_before_invasion <- 4
## number of generations after second invasion (if applicable)
generations_after_invasion <- number_generations - generations_before_invasion

## SIMULATION PARAMETERS ##
## number of monte carlo replicates
num_trials <- 100
## number of cores
num_cores <- 1
## set up parallel backend
registerDoMC(cores = num_cores)

## get name of the simulation (worker reproduction and founder status)
reproduction_mode <- if (worker_reproduction_status) "worker_reproduction" else "no_worker_reproduction"
founder_mode <- if (multiple_founder_status) "multiple_founder" else "single_founder"
sim_name <- paste(reproduction_mode, founder_mode, sep = "_")

## get and create directories, load functions file
working_directory <- getwd()
bash_date <- format(Sys.time(),format="%Y%m%d")
function_file <- sprintf("%s/QL_reproduction_functions.R", working_directory)
dir.create(file.path(working_directory, "/user_data/"), showWarnings = FALSE)
output_directory <- sprintf("%s/user_data/", getwd())
source(function_file)

results <- foreach (loop = 1:num_trials, .combine = rbind) %dopar% { 

    repeat { ## if the population goes extinct, I want to repeat the simulation (to standarize the number of "successful" invasions  

        temp_list <- list()

        ## initialise objects to store information about the simulation
        queen_allele_frequencies <- matrix(nrow = 0, ncol = number_alleles)
        worker_laid_allele_distribution <- matrix(nrow = 0, ncol = number_alleles)
        queen_laid_allele_distribution <- matrix(nrow = 0, ncol = number_alleles)
        simulation_extinction_status <- NULL ## 0 for extinct; 1 for no extinction
        generation_all_QR_colonies_lost <- NULL
        population_size <- NULL

        ## initialise invading population
        population <- initialiseInvadingColony(number_alleles, initial_distribution_alleles, number_drone_matings, cost_homozygosity)

        ## INVASION ##
        for (generation_loop in 1:number_generations){

            ## iterate one generation
            list_output <- generateNewPopulation(population, number_alleles, number_drone_matings, average_swarms,
                                  prob_queen_survives, cost_homozygosity, QL_drone_production,
                                  worker_reproduction_status)
            ## reform population
            population <- matrix(unlist(list_output[1]), nrow = length(unlist(list_output[1])) / (number_alleles + 4), ncol = number_alleles + 4)
            
            ## test for population extinction
            if ( isColonyExtinct(population, number_alleles) ){ 
                break ## breaks out of for 1:number_generations loop
            }

            ## record simulation data (queen alleles, worker-laid drone alleles, queen-laid drone alleles, population size)
            ## normalised queen allele frequencies
            tabulated_queen_alleles <- tabulate(bin = population[population[ , number_alleles + 4] == 1, 1:2], nbins = number_alleles)
            normalised_queen_alleles <-  tabulated_queen_alleles / sum(tabulated_queen_alleles)
            
            queen_allele_frequencies <- rbind( queen_allele_frequencies, normalised_queen_alleles)
            ## check to make sure queen_allele_frequencies are properly normalised
            stopifnot(abs(1 - sum(normalised_queen_alleles)) < 0.0000001)
            ## record drone allele distributions
            queen_laid_allele_distribution <- rbind( queen_laid_allele_distribution, unlist(list_output[2]) ) ## queen-laid drones
            worker_laid_allele_distribution <- rbind( worker_laid_allele_distribution, unlist(list_output[3]) ) ## worker-laid drones
            ## record number of colonies left
            population_size[length(population_size) + 1] <- getNumberOfColonies(population)
            
            ## condition for a new colony to invade from source population
            ## (this is the only difference between the single founder and multiple founder scenarios)
            if ( multiple_founder_status && generation_loop == generations_before_invasion ){
                ## add new colony from source population to invasive population
                population <- rbind(population,
                                    initialiseInvadingColony(number_alleles, initial_distribution_alleles, number_drone_matings, cost_homozygosity) )
            }
            
        }
        
        ## population is either extinct, or viable but finished number_generations
        if ( isColonyExtinct(population, number_alleles) ) { ## population is extinct
            ## record number of colonies left
            population_size[length(population_size) + 1] <- getNumberOfColonies(population)
            ## record number of generations before extinction
            generation_all_QR_colonies_lost <- generation_loop
            ## mark simulation status as "extinct"
            simulation_extinction_status <- 0
            ## simulation will "repeat"
        } else { ## population finished number_generations (don't need to record population size (done above)
            generation_all_QR_colonies_lost <- generation_loop ## will record number_generations
            ## mark simulation status as "not extinct"
            simulation_extinction_status <- 1
            break ## transfer control to foreach loop (don't want to repeat)
        }
        
    } ## end of the repeat loop
    
    ## record simulation 
    temp_list[[1]] <- queen_allele_frequencies
    temp_list[[2]] <- queen_laid_allele_distribution
    temp_list[[3]] <- worker_laid_allele_distribution
    temp_list[[4]] <- population_size
    temp_list[[5]] <- simulation_extinction_status
    temp_list[[6]] <- generation_all_QR_colonies_lost
    
    results <- temp_list
    
} ## end of the foreach loop

filename <- sprintf("%s/%s_%s_na%d_ch%d_as%d_pq%d_dm%d_qdp%d.RData",
                    output_directory,bash_date,sim_name,number_alleles,cost_homozygosity*100,average_swarms,probability_QL_colony*100,number_drone_matings, QL_drone_production)
save(results, file = filename)


