#clear all
rm(list = ls())

#simulation name
sim_name <- "QL_no_reproduction_multiple_founder"

# get and create directories, load functions file
working_directory <- getwd()
bash_date <- format(Sys.time(),format="%Y%m%d")
function_file <- sprintf("%s/%s_%s_functions.R", working_directory, bash_date, sim_name)
dir.create(file.path(working_directory, "/user_data/"), showWarnings = FALSE)
output_directory <- sprintf("%s/user_data/", getwd())
source(function_file)

#LOAD PACKAGES
require("doMC")
require("foreach")

# PARAMETERS
#starting number of colonies
N_starting_population <- 1 
#number of alleles in population
number_alleles <- 15
#number of drones each queen mates with
number_drone_matings <- 25
#cost to colony of producing diploid (dead) males
cost_homozoygosity <- 0.05
#number of generations before termination
generations_before_invasion <- 4
## total number of generations
number_generations <- 6
generations_after_invasion <- number_generations - generations_before_invasion
#proportion of swarms that become queen-less
probability_QL_colony <- 0.3
#probability that queen AND colony survives each generation
prob_queen_survives <- 1 - probability_QL_colony
#mean swarms per colony
average_swarms <- 4
#counter to record generation
counter <- 1
#initial distribution alleles (assume alleles are at balancing selection equilibrium - i.e. even frequency)
initial_distribution_alleles <- rep(1/number_alleles,number_alleles) 
#number of monte carlo replicates
num_trials <- 100
#number of cores
num_cores <- 2
#set seed (doRNG is reproducible)
seed <- 1
#set up parallel backend
registerDoMC(cores = num_cores)
set.seed(seed)

results <- foreach (loop = 1:num_trials, .combine = rbind) %dopar% { 
  extinction_counter = 0
  repeat { #the population can go extinct - if this occurs, I want to repeat the simulation 
    #initialise objects to store information about the simulation
    temp_list <- list()
    queen_allele_frequencies <- matrix(numeric((number_generations + 1) * number_alleles),
                                       nrow = number_generations + 1,ncol = number_alleles)
    drone_allele_frequencies <- matrix(numeric((number_generations + 1) * number_alleles),
                                       nrow = number_generations + 1,ncol = number_alleles)
    original_queen_survival <- numeric(number_generations + 1)
    
    population_size <- NULL
    
    N <- N_starting_population
    
    counter <- 1
    
    #produce population matrix 
    population <- matrix(numeric(N * (2 + number_alleles + 2)),nrow = N, ncol = 2 + number_alleles + 2)
    
    #choose queen alleles (cols 1 and 2 of population matrix)
    for (ii in 1:N){
      population[ii,1:2] <- sample(1:number_alleles, 2, replace = TRUE, prob = initial_distribution_alleles)
      
      while (population[ii,1] == population[ii,2]){ #if homozygous, 'kill' and choose again
        population[ii,1:2] <- sample(1:number_alleles, 2, replace = TRUE, prob = initial_distribution_alleles)
      }
    }
    
    #choose which drones each queen mates with 
    sampled_drones <- 
      matrix(sample(1:number_alleles,N * number_drone_matings,replace = TRUE, 
                    prob = initial_distribution_alleles), nrow = N, ncol = number_drone_matings)
    
    #transform these alleles into proportions
    drone_proportions <- matrix(numeric(N * number_alleles),nrow = N, ncol = number_alleles)
    for (ii in 1:N){
      drone_proportions[ii,] <- tabulate(bin = sampled_drones[ii,],nbins = number_alleles) / number_drone_matings
    }
    
    #cols 3:(number_alleles + 2) represent the spermathecal contents
    population[,3:(number_alleles + 2)] <- drone_proportions
    
    #calculate colony fitnesses (for queen-less colonies, fitness only affects drone production)
    for (ii in 1:N){
      # determine proportion of homozygosity by checking each queen allele and multiplying it by the corresponding drone allele
      # each queen allele contributes half the total homozygosity and total homozygosity is multipled by its cost
      population[ii,number_alleles + 3] <- 1 - (((0.5 * population[ii,population[ii,1] + 2]) + 
                                                   (0.5 * population[ii,population[ii,2] + 2])) * cost_homozoygosity)
      
      population[ii,number_alleles + 4] <- 1 #first colony is queenright (1 = QR, 0 = QL, delete row = dead)
    }
    
    #record information about the simulation
    queen_allele_frequencies[counter,] <- sort(tabulate(bin = population[,1:2],nbins = number_alleles) /
                                                 (N * 2), decreasing = TRUE) #counts, ordered from highest freq to lowest
    
    #drone alleles are stored as proportions not counts so use colMeans
    #if N == 1, colMeans throws an error (array becomes vector and it requires 2 dim object) so use drop = FALSE
    drone_allele_frequencies[counter,] <- sort(colMeans(population[,3:(number_alleles + 2), drop = FALSE]),decreasing = TRUE) 
    
    population_size[counter] <- N
    
    original_colony <- population[1,]
    
    original_queen_survival[counter] <- identical(population[1,],original_colony)
    
    counter <- counter + 1
    
    for (generationloop in 1:generations_before_invasion){
      
      #logistic growth rate (depends on population size, we include both QL and QR colonies)
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
        
        break #breaks out of for 1:generations_before_invasion loop
        
      } else { #not extinct, proceed normally
        
        #no longer need the old population size 
        N <- nrow(population)
        
        #record information about the simulation
        queen_allele_frequencies[counter,] <- sort(tabulate(bin = population[,1:2],nbins = number_alleles) /
                                                     (N * 2), decreasing = TRUE) #counts, ordered from highest freq to lowest
        
        #drone alleles are stored as proportions not counts so use colMeans
        #if N == 1, colMeans throws an error (array becomes vector and it requires 2 dim object) so use drop = FALSE
        drone_allele_frequencies[counter,] <- sort(colMeans(population[,3:(number_alleles + 2), drop = FALSE]),decreasing = TRUE) 
        
        population_size[counter] <- N
        
        original_queen_survival[counter] <- identical(population[1,],original_colony)
        
        counter <- counter + 1
      } 
    }
    
    ###### NEW COLONY INVADES ######
    
    if (population != 'extinct'){
      #add row of zeros to population matrix
      population <- rbind(population,numeric(2 + number_alleles + 1))
      #adjust population size (now N + 1 not N)
      N <- N + 1
      
      #choose queen alleles (cols 1 and 2 of population matrix)
      population[N,1:2] <- sample(1:number_alleles, 2, replace = TRUE, prob = initial_distribution_alleles)
      
      while (population[N,1] == population[N,2]){ #if homozygous, 'kill' and choose again
        population[N,1:2] <- sample(1:number_alleles, 2, replace = TRUE, prob = initial_distribution_alleles)
      }
      
      #choose which drones each queen mates with 
      sampled_drones <- sample(1:number_alleles,number_drone_matings,replace = TRUE, prob = initial_distribution_alleles)
      
      #transform these alleles into proportions
      drone_proportions<- tabulate(bin = sampled_drones,nbins = number_alleles) / number_drone_matings
      
      #cols 3:(number_alleles + 2) represent the spermathecal contents
      population[N,3:(number_alleles + 2)] <- drone_proportions
      
      #calculate invasive colony fitness
      population[N,number_alleles + 3] <- 1 - (((0.5 * population[N,population[N,1] + 2]) + 
                                                  (0.5 * population[N,population[N,2] + 2])) * cost_homozoygosity)
      
    }
      
    for (generationloop in 1:generations_after_invasion){
      
      #is the population extinct?
      if (population == 'extinct'){ 
        
        break #breaks out of for 1:number_generations loop
        
      }
      
      #logistic growth rate (depends on population size, we include both QL and QR colonies)
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
        
        break #breaks out of for 1:generations_before_invasion loop
        
      } else { #not extinct, proceed normally
        
        #no longer need the old population size 
        N <- nrow(population)
        
        #record information about the simulation
        queen_allele_frequencies[counter,] <- sort(tabulate(bin = population[,1:2],nbins = number_alleles) /
                                                     (N * 2), decreasing = TRUE) #counts, ordered from highest freq to lowest
        
        #drone alleles are stored as proportions not counts so use colMeans
        #if N == 1, colMeans throws an error (array becomes vector and it requires 2 dim object) so use drop = FALSE
        drone_allele_frequencies[counter,] <- sort(colMeans(population[,3:(number_alleles + 2), drop = FALSE]),decreasing = TRUE) 
        
        population_size[counter] <- N
        
        original_queen_survival[counter] <- identical(population[1,],original_colony)
        
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
    
  } #this marks the end of the repeat loop
  
  #record simulation 
  
  temp_list[[1]] <- queen_allele_frequencies
  temp_list[[2]] <- drone_allele_frequencies
  temp_list[[3]] <- population_size
  temp_list[[4]] <- extinction_counter
  temp_list[[5]] <- sum(original_queen_survival)
  
  results <- temp_list
  
} #this marks the end of the foreach loop

filename <- sprintf("%s/%s_%s_na%d_ch%d_as%d_pq%d_dm%d.RData",
                    output_directory,bash_date,sim_name,number_alleles,cost_homozoygosity*100,average_swarms,probability_QL_colony*100,number_drone_matings)
save(results, file = filename)

