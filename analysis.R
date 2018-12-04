## Joshua Christie (joshua.christie@sydney.edu.au)

rm(list=ls())

## function to get parameter values from script name
getParameterValue <- function(parameter, string){
    regex_string <- sprintf('^.*%s\\s*|\\s*(_|\\.).*$', parameter)
    return(as.double(gsub(regex_string, '', string)))
}

## parameter values
num_cols <- 49
num_gen <- 6
counter <- 1
dm <- 25
## script is set to analyse the simulation data in Gloag et al. Molecular Ecology, which are found in manuscript_data/
path_to_folder <- sprintf('%s/manuscript_data/', getwd())
## uncomment the line below if you instead want to analyse your own simulations in user_data/
## path_to_folder <- sprintf('%s/user_data/', getwd())

file_names <- list.files(path=path_to_folder,pattern=".*.RData")

summary <- data.frame()

for (filename in file_names){
    
    load(sprintf('%s%s', path_to_folder, filename))
    
    summary <- rbind(summary, numeric(num_cols))
    
    summary[counter, 1] <- sub("_na.*", "", substring(filename,10)) ## simulation type (worker reproduction or not; multiple founder or not)
    summary[counter, 2] <- getParameterValue('ch', filename) / 100 ## cost of homozygosity
    summary[counter, 3] <- getParameterValue('na', filename) ## number of alleles in the source population
    summary[counter, 4] <- getParameterValue('as', filename) ## average swarms left by a colony with relative fitness of 1
    summary[counter, 5] <- getParameterValue('pq', filename) / 100 ## probability of queenlessness
    ## competitiveness of worker-laid drones relative to queen-laid drones (only relevant for simulations with worker reproduction)
    summary[counter, 6] <- if (summary[counter, 1] == "no_worker_reproduction_single_founder" || summary[counter, 1] == "no_worker_reproduction_multiple_founder") "NA" else getParameterValue('qdp', filename) / 10

    ## data in loaded file are organised as follows:
    
    ## results[[1]] <- queen_allele_frequencies
    ## results[[2]] <- queen_laid_allele_distribution
    ## results[[3]] <- worker_laid_allele_distribution
    ## results[[4]] <- population_size
    ## results[[5]] <- simulation_extinction_status
    ## results[[6]] <- generation_all_QR_colonies_lost
    
    ## set up some arrays that will store various metrics for each simulation (rows) in each generation (cols)
    queen_alleles_matrix <- matrix(nrow=0,ncol=num_gen)
    worker_laid_alleles_matrix <- matrix(nrow=0,ncol=num_gen)
    queen_laid_alleles_matrix <- matrix(nrow=0,ncol=num_gen)
    DMP_combined <- matrix(nrow=0,ncol=num_gen)
    difference_queen_worker_drone <- matrix(nrow=0,ncol=num_gen)
    mean_freq_difference_queen_worker_drone <- {} ## record mean frequency of drone alleles in worker-laid but not queen-laid fraction
    population_size_vector <- {} ## this will record population size for each simulation but for generation 6 only
    
    ## get indices of those simulations that completed without the population going extinct
    sims_with_no_extinction <- as.logical(unlist(results[,5]))
    ## get number of simulations
    num_sims <- length(results[,1]) ## 5000 for the data in the manuscript
    ## loop through each simulation
    for (i in 1:num_sims){
        
        ## record data for sims in which the population didn't go extinct (I analyse extinction data separately below)
        if (sims_with_no_extinction[i]){
            
            ## alleles carried by queens and the two categories of drones
            queen_alleles <- results[[i,1]]
            queen_laid_alleles <- results[[i,2]] ## queen-laid drones
            worker_laid_alleles <- results[[i,3]] ## worker-laid drones
            
            ## get population size at gen 6
            population_size_vector[length(population_size_vector) + 1] <- results[[i,4]][num_gen]

            ## set up some temporary vectors (these will be rbinded to the arrays above)
            diff_temp <- numeric(num_gen)
            DMP_combined_temp <- numeric(num_gen)
            queen_alleles_temp <- numeric(num_gen)
            worker_laid_alleles_temp <- numeric(num_gen)
            queen_laid_alleles_temp <- numeric(num_gen)

            ## now loop through each generation
            for (j in 1:num_gen){

                ## get number of alleles present for queens and workers
                queen_alleles_temp[j] <- sum(queen_alleles[j,]>0)
                queen_laid_alleles_temp[j] <- sum(queen_laid_alleles[j,]>0)
                worker_laid_alleles_temp[j] <- sum(worker_laid_alleles[j,]>0)
                
                ## get indices of queen-laid drone alleles present at non-zero frequency
                queen_laid_logical <- queen_laid_alleles[j, ] > 0
                worker_laid_logical <- worker_laid_alleles[j, ] > 0 ## ditto for worker-laid drones
                ## get indices of alleles present in worker-laid drones but not in queen-laid drones
                difference_logical <- !queen_laid_logical & worker_laid_logical
                ## get number of alleles fitting this criterion
                diff_temp[j] <- sum(difference_logical)
                ## get mean frequency of alleles fitting this criterion
                if (j == num_gen ){ ## only want gen 6
                    mean_freq_difference_queen_worker_drone[length(mean_freq_difference_queen_worker_drone) + 1] <-
                        mean(worker_laid_alleles[j, difference_logical])
                }
                
                stopifnot(abs(sum(queen_laid_alleles[j, ] + worker_laid_alleles[j, ])- 1) < 0.00000001)
                
                ## get DMP level
                DMP_combined_temp[j] <- sum(queen_alleles[j, ] * (queen_laid_alleles[j, ] + worker_laid_alleles[j, ]))

            }
            
            ## rbind temporary vectors to relevant array
            difference_queen_worker_drone <- rbind(difference_queen_worker_drone, diff_temp)
            DMP_combined <- rbind(DMP_combined, DMP_combined_temp)
            queen_alleles_matrix <- rbind(queen_alleles_matrix, queen_alleles_temp)
            queen_laid_alleles_matrix <- rbind(queen_laid_alleles_matrix, queen_laid_alleles_temp)
            worker_laid_alleles_matrix <- rbind(worker_laid_alleles_matrix, worker_laid_alleles_temp)
        }
        
    }
    ## get mean queen allele diversity among the simulations for each generation
    summary[counter,7] <- mean(queen_alleles_matrix[,1])
    summary[counter,8] <- mean(queen_alleles_matrix[,2])
    summary[counter,9] <- mean(queen_alleles_matrix[,3])
    summary[counter,10] <- mean(queen_alleles_matrix[,4])
    summary[counter,11] <- mean(queen_alleles_matrix[,5])
    summary[counter,12] <- mean(queen_alleles_matrix[,6])
    ## get sd and median for the final gen
    summary[counter,13] <- sd(queen_alleles_matrix[,6])
    summary[counter,14] <- median(queen_alleles_matrix[,6])
    ## get mean diversity of queen-laid drones
    summary[counter,15] <- mean(queen_laid_alleles_matrix[,1])
    summary[counter,16] <- mean(queen_laid_alleles_matrix[,2])
    summary[counter,17] <- mean(queen_laid_alleles_matrix[,3])
    summary[counter,18] <- mean(queen_laid_alleles_matrix[,4])
    summary[counter,19] <- mean(queen_laid_alleles_matrix[,5])
    summary[counter,20] <- mean(queen_laid_alleles_matrix[,6])
    ## sd and median for final gen
    summary[counter,21] <- sd(queen_laid_alleles_matrix[,6])
    summary[counter,22] <- median(queen_laid_alleles_matrix[,6])
    ## get mean diversity of worker-laid drones
    summary[counter,23] <- mean(worker_laid_alleles_matrix[,1])
    summary[counter,24] <- mean(worker_laid_alleles_matrix[,2])
    summary[counter,25] <- mean(worker_laid_alleles_matrix[,3])
    summary[counter,26] <- mean(worker_laid_alleles_matrix[,4])
    summary[counter,27] <- mean(worker_laid_alleles_matrix[,5])
    summary[counter,28] <- mean(worker_laid_alleles_matrix[,6])
    ## sd and median for final gen
    summary[counter,29] <- sd(worker_laid_alleles_matrix[,6])
    summary[counter,30] <- median(worker_laid_alleles_matrix[,6])
    ## get mean difference in allele diversity between queen-laid and worker-laid drones
    summary[counter,31] <- mean(difference_queen_worker_drone[,1])
    summary[counter,32] <- mean(difference_queen_worker_drone[,2])
    summary[counter,33] <- mean(difference_queen_worker_drone[,3])
    summary[counter,34] <- mean(difference_queen_worker_drone[,4])
    summary[counter,35] <- mean(difference_queen_worker_drone[,5])
    summary[counter,36] <- mean(difference_queen_worker_drone[,6])
    ## get mean allele freq in fraction of worker-laid alleles that don't appear in the queen-laid fraction
    summary[counter,37] <- mean(mean_freq_difference_queen_worker_drone, na.rm = TRUE)
    ## get proportion of extinctions
    summary[counter,38] <- 1-sum(unlist(results[,5]))/num_sims

    ## mean generation in which extinctions occur
    gen_QR_lost <- unlist(results[,6][!sims_with_no_extinction])
    ## get mean generation in which extinctions occur
    summary[counter,39] <- mean(gen_QR_lost)
    ## get proportion of extinctions that occur in generation 1
    summary[counter,40] <- sum(gen_QR_lost == 1)/length(gen_QR_lost)
    ## get mean population size at gen 6
    summary[counter,41] <- mean(population_size_vector)
    ## get mean DMP production for each generation
    summary[counter,42] <- mean(DMP_combined[ , 1])
    summary[counter,43] <- mean(DMP_combined[ , 2])
    summary[counter,44] <- mean(DMP_combined[ , 3])
    summary[counter,45] <- mean(DMP_combined[ , 4])
    summary[counter,46] <- mean(DMP_combined[ , 5])
    summary[counter,47] <- mean(DMP_combined[ , 6])
    ## sd and range for gen 6
    summary[counter,48] <- sd(DMP_combined[ , 6])
    summary[counter,49] <- sprintf('%1.2f--%1.2f', range(DMP_combined[ , 6])[1],range(DMP_combined[ , 6])[2])

    counter <- counter + 1
    
}

colnames(summary) <- c("sim type", "cost homozygosity",	"alleles", "swarm", "QLness rate", "QL drone competitiveness",  "mean queen alleles @ gen 1","mean queen alleles @ gen 2","mean queen alleles @ gen 3","mean queen alleles @ gen 4","mean queen alleles @ gen 5", "mean queen alleles @ gen 6", "sd @ gen 6", "median @ gen 6", "mean queen-laid drone alleles @ gen 1", "mean queen-laid drone alleles @ gen 2", "mean queen-laid drone alleles @ gen 3", "mean queen-laid drone alleles @ gen 4", "mean queen-laid drone alleles @ gen 5", "mean queen-laid drone alleles @ gen 6", "sd @ gen 6", "median @ gen 6", "mean worker-laid drone alleles @ gen 1", "mean worker-laid drone alleles @ gen 2", "mean worker-laid drone alleles @ gen 3", "mean worker-laid drone alleles @ gen 4", "mean worker-laid drone alleles @ gen 5", "mean worker-laid drone alleles @ gen 6", "sd @ gen 6", "median @ gen 6", "mean alleles carried in worker-laid fraction not found in queen-laid fraction @ gen 1", "mean alleles carried in worker-laid fraction not found in queen-laid fraction @ gen 2", "mean alleles carried in worker-laid fraction not found in queen-laid fraction @ gen 3","mean alleles carried in worker-laid fraction not found in queen-laid fraction @ gen 4","mean alleles carried in worker-laid fraction not found in queen-laid fraction @ gen 5","mean alleles carried in worker-laid fraction not found in queen-laid fraction @ gen 6", "mean frequency of alleles carried in worker-laid fraction not found in queen-laid fraction @ gen 6","proportion of extinctions", "mean gen in which extinction occurs (all QR colonies lost)", "proportion of extinctions that were @ gen 1", "population size @ gen6", "DMP @ gen 1", "DMP @ gen 2", "DMP @ gen 3", "DMP @ gen 4", "DMP @ gen 5", "DMP at gen 6", "sd @ gen6", "range @ gen6")

output_filename <- sprintf('%sdata_summary.csv', path_to_folder)

write.csv(summary, output_filename, row.names=FALSE)
