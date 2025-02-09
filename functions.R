initialiseInvadingColony <- function(number_alleles, initial_distribution_alleles, number_drone_matings, cost_homozygosity){
    ## set up invading colony with queen alleles and drone matings from a source population    

    ## initialise new colony
    invading_colony <- matrix(numeric(2 + number_alleles + 2),nrow = 1, ncol = 2 + number_alleles + 2)

    ## choose queen alleles of invader from the source population
    invading_colony <- initialiseQueenAllelesFromSourcePop(invading_colony, number_alleles, initial_distribution_alleles)

    ## choose which drones the invading queen mates with (also from source population, since she arrives already mated)
    invading_colony <- chooseDroneAlleles(invading_colony, colony_ID = 1, number_alleles, number_drone_matings, initial_distribution_alleles)
                                              
    ## calculate fitness of invading colony
    invading_colony <- calculateColonyFitness(invading_colony, colony_ID = 1, number_alleles, cost_homozygosity)
    
    ## assume colony remains queenright
    invading_colony <- setColonyQueenStatus(invading_colony, number_alleles, colony_ID = 1, queen_status = 1)

    return(invading_colony)
}

initialiseQueenAllelesFromSourcePop <- function(population, number_alleles, initial_distribution_alleles){
    ## choose alleles of invading queen from the source population (cols---1:2 of population)

    repeat{ ## females must be heterozygous
        ## choose queen alleles from source population
        population[1, 1:2] <- sample(1:number_alleles, 2, replace = TRUE, prob = initial_distribution_alleles)
        if ( population[1, 1] != population[1, 2] ){ ## heterozygous
            break
        }
        
    }
    return(population)
}

chooseDroneAlleles <- function(population, colony_ID, number_alleles, number_drone_matings, allele_distribution){
    ## choose drones that mate with the queen, turn allele counts into proportions, and add the queen's
    ## spermathecal contents to the colony_IDth row of the population (cols---3:(number_alleles + 2))
    
    ## sample the number_alleles in the population, in proportion to the drone allele_distribution, number_drone_matings times
    sampled_drones <- sample(1:number_alleles, number_drone_matings, replace = TRUE, prob = allele_distribution)
    ## transform these alleles into proportions
    drone_proportions <- tabulate(bin = sampled_drones, nbins = number_alleles) / number_drone_matings
    ## cols 3:(number_alleles + 2) represent the spermathecal contents
    population[colony_ID, 3:(number_alleles + 2)] <- drone_proportions
    
    return(population)
}

calculateColonyFitness <- function(population, colony_ID, number_alleles, cost_homozygosity){
    ## add the colony fitness of the colony_IDth row of the population (col---number_alleles + 3)

    ## for QR colonies, fitness affects production of drones and daughter queens
    ## for QL colonies, fitness only affects drone production (and only if worker_reproduction_status is TRUE)
    queen_allele_1_ID <- population[colony_ID, 1]
    queen_allele_2_ID <- population[colony_ID, 2]
    ## determine proportion of homozygosity by multiplying each queen allele frequency (0.5) with the corresponding drone allele
    ## corresponding drone alleles are shifted 2 cols to the right (since the first 2 cols store the queen's genotype)
    homozygosity_level <- 0.5 * population[colony_ID, queen_allele_1_ID + 2] + 0.5 * population[colony_ID, queen_allele_2_ID + 2]
    homozygosity_fitness_cost <- homozygosity_level * cost_homozygosity
    ## set colony fitness
    population[colony_ID, number_alleles + 3] <- 1 - homozygosity_fitness_cost

    return(population)
}

setColonyQueenStatus <- function(population, number_alleles, colony_ID, queen_status){
    ## set colony status (1 = QR, 0 = QL) of colony_IDth row of population (col---number_alleles + 4)

    population[colony_ID, number_alleles + 4] <- queen_status 
    return(population)        
}

generateNewPopulation <- function(population, number_alleles, number_drone_matings, average_swarms,
                                  prob_queen_survives, cost_homozygosity, QL_drone_production,
                                  worker_reproduction_status){
    ## calls functions to iterate the population through a full generation
    
    ## initialise matrix for new population and vectors for drone alleles
    new_population <- matrix(nrow = 0, ncol = 2 + number_alleles + 2)
    worker_laid_drone_alleles <- numeric(number_alleles)
    queen_laid_drone_alleles <- numeric(number_alleles)
    
    for ( i in 1:getNumberOfColonies(population) ){ ## loop through current colonies
        
        if ( isColonyQR(colony_ID = i, number_alleles, population) ){ ## is colony QR?
            ## the queen of a QR colony does not remate, so copy colony i exactly
            new_population <- rbind( new_population, population[i, ] )
            colony_index <- getNumberOfColonies(new_population) ## will now operate on new_population[colony_index]
            ## now the only consideration left is whether the colony remains QR
            ## if not, we need to change its status to QL
            colony_i_status <- as.integer(
                testAgainstRandomNumber( prob_queen_survives * getColonyFitness(colony_index, number_alleles, new_population) ) )
            new_population <- setColonyQueenStatus(new_population, number_alleles, colony_index, colony_i_status)
            
            ## QR colony produces drones (no difference between those remaining QR and those becoming QL--we assume the queen is lost after drone production)
            queen_laid_drone_alleles <- produceDronesQueenright(new_population, number_alleles, colony_index, queen_laid_drone_alleles)
            ## QR colonies can produce daughter queens (even if they become QL because queen is lost after daughter queens have been raised)
            number_daughter_colonies <- rpois( n = 1, lambda = average_swarms * getColonyFitness(colony_index, number_alleles, new_population) )
            
            if ( as.logical(number_daughter_colonies) ){ ## at least 1 daughter colony is produced
                ## loop through these daughter colonies
                for (j in 1:number_daughter_colonies){
                    ## set up daughter colony
                    new_population <- rbind( new_population, setupDaughterColony(new_population, colony_index, number_alleles, prob_queen_survives) )
                }
            }

        } else if ( !isColonyQR(colony_ID = i, number_alleles, population) ) { ## it is QL (and became so last generation)
            ## QL colony does not survive (so don't add to new_population) nor does it produce daughter queens
            ## if workers can reproduce, colony produces drones via workers (proportional to QL_drone_production)
            if ( worker_reproduction_status ){ ## workers can reproduce
                worker_laid_drone_alleles <- produceDronesQueenless(population, number_alleles, colony_ID = i, worker_laid_drone_alleles, QL_drone_production)
            }
            
        } else {
            ## shouldn't execute
            stopifnot(FALSE)
        }

    }

    ## last step is for daughter queens that survive to mate (giving spermathecal contents), and from this, set colony fitness
    ## first, normalise the worker- and queen-laid drones so that they sum to 1
    drone_allele_sum <- sum(worker_laid_drone_alleles + queen_laid_drone_alleles)
    worker_laid_drone_alleles <- worker_laid_drone_alleles / drone_allele_sum
    queen_laid_drone_alleles <- queen_laid_drone_alleles / drone_allele_sum
    
    ## loop through each colony in new_population, find surviving daughter queens (who will mate), and set spermathecal contents and colony fitness
    for ( i in 1:getNumberOfColonies(new_population) ){

        if ( sum( new_population[i, 3:(number_alleles + 3)] ) == 0 ){ ## surviving daughter queens
            new_population <- chooseDroneAlleles(new_population, colony_ID = i, number_alleles, number_drone_matings, worker_laid_drone_alleles + queen_laid_drone_alleles)
            new_population <- calculateColonyFitness(new_population, colony_ID = i, number_alleles, cost_homozygosity)
        }

    }

    ## sanity checks
    ## workers should not reproduce if worker reproduction is turned off
    if (!worker_reproduction_status) stopifnot( !sum(worker_laid_drone_alleles) )
    ## drone alleles should be normalised
    stopifnot( abs(1 - sum(worker_laid_drone_alleles + queen_laid_drone_alleles)) < 0.0000001) 
    
    ## in addition to population, I also want to return the drone distributions (for analysis)
    list_output <- list(new_population, queen_laid_drone_alleles, worker_laid_drone_alleles)
    return(list_output)
    
}

setupDaughterColony <- function(population, mother_of_daughter_ID, number_alleles, prob_queen_survives){
    ## get vector that I will rbind to new_population
    daughter_queen <- numeric(2 + number_alleles + 2)
    ## check whether the daughter queen dies or suvives
    if ( testAgainstRandomNumber(prob_queen_survives) ){ ## daughter queen survives
        ## get distribution of alleles in her mother's spermatheca
        spermathecal_contents <- population[mother_of_daughter_ID, 3:(number_alleles + 2)]
        ## denote daughter queen's colony as QR (since she survives)
        daughter_queen[number_alleles + 4] <- 1
        
        repeat { ## females must be heterozygous
            ## allele of daughter queen that comes from her mother's genotype
            daughter_queen[1] <- population[mother_of_daughter_ID, sample(1:2, 1)]
            ## allele from daughter queen's mother's spermatheca
            daughter_queen[2] <- sample(1:number_alleles, 1, prob = spermathecal_contents)
            ## break out once daughter queen is heterozygous
            if ( daughter_queen[1] != daughter_queen[2] ){
                break
            }
        }
        
    } else {
        ## daughter queen dies, colony becomes QL---workers for this colony come from the daughter queen's mother
        ## to capture this, make the "queen" of this QL colony the daughter queen's mother
        ## need queen genotype (1:2), spermatheca (3:(number_alleles + 2), and colony fitness (number_alleles + 3)
        ## QL/QR indicator (number_alleles + 4) is initialised as zero (indicating QL), so leave this alone
        daughter_queen[1:(number_alleles + 3)] <- population[mother_of_daughter_ID, 1:(number_alleles + 3)]
    }
    
    return(daughter_queen)
    ## when the daughter colony remains QR, we still need spermathecal contents
    ## and colony fitness (for these, we must wait until daughter's mating flight)
    ## but when the daughter colony becomes immediately QL, the entry is complete
}

produceDronesQueenright <- function(population, number_alleles, colony_ID, queen_laid_drone_alleles){
    ## QR colonies produce drones (proportional to colony fitness). Drones come from queen's two alleles
    ## index into the queen alleles, adding 0.5 (frequency of each allele) multiplied by colony fitness
    
    ## initialise vector to store drone alleles from colony i 
    colony_i_drone_alleles <- numeric(number_alleles)
    ## add alleles from queen's genotype (adjusted by colony fitness)
    colony_i_drone_alleles[ population[colony_ID, 1:2] ] <- getColonyFitness(colony_ID, number_alleles, population) * 0.5
    ## add to queen_laid_drone_alleles
    queen_laid_drone_alleles <- queen_laid_drone_alleles + colony_i_drone_alleles

    return(queen_laid_drone_alleles)
}

produceDronesQueenless <- function(population, number_alleles, colony_ID, worker_laid_drone_alleles, QL_drone_production){
    ## QL colonies produce drones (proportional to colony fitness AND QL_drone_production)
    ## drones come from workers, which means 0.5 from queen's two alleles (0.25 each) and 0.5 from queen's spermatheca
    ## the latter is the contribution of the fathers of the workers
    
    ## initialise vector to store drone alleles from colony i 
    colony_i_drone_alleles <- numeric(number_alleles)
    ## add alleles from queen's genotype (adjusted by colony fitness and multiplied by 0.25 instead of 0.5)
    colony_i_drone_alleles[ population[colony_ID, 1:2] ] <- getColonyFitness(colony_ID, number_alleles, population) * 0.25
    ## get spermathecal contents of queen
    spermathecal_contents <- population[colony_ID, 3:(number_alleles + 2)]
    ## add alleles from queen's spermatheca (contribution of workers' fathers)
    colony_i_drone_alleles <- colony_i_drone_alleles + spermathecal_contents *
        getColonyFitness(colony_ID, number_alleles, population) * 0.5
    ## account for different output of drones from QL colonies (compared to QR colonies)
    colony_i_drone_alleles <- colony_i_drone_alleles * QL_drone_production
    ## add to worker_laid_drone_alleles
    worker_laid_drone_alleles <- worker_laid_drone_alleles + colony_i_drone_alleles

    return(worker_laid_drone_alleles)
}

isColonyQR <- function(colony_ID, number_alleles, population){
    
    return( as.logical(population[colony_ID, number_alleles + 4]) )
    
}

testAgainstRandomNumber <- function(threshold){

    return( runif(1) < threshold )

}

getColonyFitness <- function(colony_ID, number_alleles, population){

    return( population[colony_ID, number_alleles + 3] )

}

isColonyExtinct <- function(population, number_alleles){
    ## test for extinction (no colonies added or any colonies added are QL)
    return( !nrow(population) || !sum( population[ , number_alleles + 4] ) )
}

getNumberOfColonies <- function(population){
    ## returns number of colonies in the population

    ## nrow doesn't work on vectors (gives NULL) but NROW treats vector as a col vector
    if ( !is.null( nrow(population) ) ) {
        number_colonies <- nrow(population)
    } else { ## if nrow returns NULL, the matrix has collapsed to a vector and there's 1 colony
        number_colonies <- 1
    }
    return(number_colonies)
}
