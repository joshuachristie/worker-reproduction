initialiseQueenAllelesFromSourcePop <- function(population, number_alleles, initialise_distribution_alleles){

    repeat{ ## females must be heterozygous
        ## choose queen alleles from source population
        population[1, 1:2] <- sample(1:number_alleles, 2, replace = TRUE, prob = initial_distribution_alleles)
        if ( population[1, 1] != population[1, 2] ){ ## heterozygous
            break
        }
        
    }
    return(population)
}

setupDaughterColony <- function(population, new_population, old_colony_ID, number_alleles, prob_queen_survives){
    ## get vector that I will rbind to new_population
    daughter_queen <- numeric(2 + number_alleles + 2)
    ## check whether the daughter queen dies or suvives
    if ( testAgainstRandomNumber(prob_queen_survives) ){ ## daughter queen survives
        ## get distribution of alleles in her mother's spermatheca
        spermathecal_contents <- population[old_colony_ID, 3:(number_alleles + 2)]
        ## denote daughter queen's colony as QR (since she survives)
        daughter_queen[number_alleles + 4] <- 1
        
        repeat { ## females must be heterozygous
            ## allele of daughter queen that comes from her mother's genotype
            daughter_queen[1] <- population[old_colony_ID, sample(1:2, 1)]
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
        daughter_queen[1:(number_alleles + 3)] <- population[old_colony_ID, 1:(number_alleles + 3)]
    }
    
    new_population <- rbind(new_population, daughter_queen)
    return(new_population)
    ## when the daughter colony remains QR, we still need spermathecal contents
    ## and colony fitness (for these, we must wait until daughter's mating flight)
    ## but when the daughter colony becomes immediately QL, the entry is complete
}

chooseDroneAlleles <- function(population, colony_ID, number_alleles, number_drone_matings, allele_distribution){
    ## choose drones that mate with the queen, turn allele IDs into proportions, and add the queen's spermathecal contents to population
    
    ## sample the number_alleles in the population, in proportion to the drone allele_distribution, number_drone_matings times
    sampled_drones <- sample(1:number_alleles, number_drone_matings, replace = TRUE, prob = allele_distribution)
    ## transform these alleles into proportions
    drone_proportions <- tabulate(bin = sampled_drones, nbins = number_alleles) / number_drone_matings
    ## cols 3:(number_alleles + 2) represent the spermathecal contents
    population[colony_ID, 3:(number_alleles + 2)] <- drone_proportions
    
    return(population)
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

calculateColonyFitness <- function(population, colony_ID, number_alleles, cost_homozygosity){
    ## for QR colonies, fitness affects production of drones and daughter queens
    ## for QL colonies, fitness only affects drone production
    queen_allele_1_ID <- population[colony_ID, 1]
    queen_allele_2_ID <- population[colony_ID, 2]
    ## determine proportion of homozygosity by multiplying each queen allele frequency (0.5) with the corresponding drone allele
    ## corresponding drone alleles are shifted 2 cols to the right (since the first 2 cols store the queen's genotype)
    homozygosity_level <- 0.5 * population[colony_ID, queen_allele_ID_1 + 2] + 0.5 * population[colony_ID, queen_allele_ID_2 + 2]
    homozygosity_cost <- homozygosity_level * cost_homozygosity
    ## set colony fitness
    population[colony_ID, number_alleles + 3] <- 1 - homozygosity_cost

    return(population)
}

setColonyQueenStatus <- function(population, number_alleles, colony_ID, queen_status){
    ## set colony status (1 = QR, 0 = QL)
    population[colony_ID, number_alleles + 4] <- queen_status 
    return(population)        
}

isColonyQR <- function(colony_ID, number_alleles, population){
    
    return( as.logical(population[colony_ID, number_alleles + 4]) )
    
}

testAgainstRandomNumber <- function(threshold){

    return( runif(1) < threshold )

}

    
                                        #determine which colonies survive intact (queen + colony)
    queenright_survival <- runif(N) < #depends on fitness of colony AND on population size
        (prob_queen_survives * population[,number_alleles + 3]) 
                                        #determine which colonies survive without queen
    queenless_survival <- rep(TRUE,N)
    
    queen_survival[queenright_survival] <- 1 #colony survives queenright
    queen_survival[!queenright_survival & !queenless_survival] <- -1 #colony dies
                                        #remainder are !queenright_survival & queenless_survival, which store 0
    
###queen_survival[ii] reads 1 when iith colony is QR, 0 when QL, -1 when dead###  
    
                                        #determine the size of the next generation (daughter swarms of QR parents ONLY + all parent QR colonies that survive (as QR or QL))
    N_new_generation <- N_daughter_colonies + sum((population[,number_alleles + 4] == 1 & queen_survival > -1))
    
                                        #at this point, the population can go extinct. If it does, I want to escape the function and
                                        # repeat the foreach loop in the main script (so that I get 500 simulations in which the population does not go extinct).
    ## ADD HERE A CONDITION THAT THERE ARE QR COLONIES (OTHERWISE POPULATION WILL GO EXTINCT THE NEXT GENERATION)
    if (N_new_generation == 0){ #population is extinct 
        
        population <- 'extinct' #for conditional statement in main script
        return(population) #exit function
        
    } else { #population is not extinct, proceed normally
        
                                        #initialise arrays
        new_population <- matrix(numeric(N_new_generation * (2 + number_alleles + 2)),
                                 nrow = N_new_generation, ncol = 2 + number_alleles + 2)
        distribution_alleles <- numeric(number_alleles)

        ## If a new colony goes QL, I need to be able to track which colony it came from
        population_original_queen_index <- numeric(N_new_generation)

        for (ii in 1:N){
                                        #is colony ii QR? If so, does it remain QR?
            if (queen_survival[ii] == 1 && population[ii,number_alleles + 4] == 1){ 
                
                new_population[swarm_counter,] <- population[ii,] #if TRUE, copy colony ii exactly

                ## store which colony this entry is derived from
                population_original_queen_index[swarm_counter] <- ii

                swarm_counter <- swarm_counter + 1 
                
                                        #QR colony produces drones
                distribution_alleles[population[ii,1:2]] <- distribution_alleles[population[ii,1:2]] + 
                    (population[ii,number_alleles + 3] * 0.5) 
                
                                        #if not, is colony ii QR but becomes QL?
            } else if (queen_survival[ii] == 0 && population[ii,number_alleles + 4] == 1){ 

                                        #copy colony ii exactly except for its QR/QL indicator column
                new_population[swarm_counter,1:(number_alleles + 3)] <- population[ii,1:(number_alleles + 3)] 
                                        #change colony indicator to QL
                new_population[swarm_counter,number_alleles + 4] <- 0

                ## store which colony this entry is derived from
                population_original_queen_index[swarm_counter] <- ii

                swarm_counter <- swarm_counter + 1 
                
                                        #QR colony produces drones
                distribution_alleles[population[ii,1:2]] <- distribution_alleles[population[ii,1:2]] + 
                    (population[ii,number_alleles + 3] * 0.5) 
                
                                        # if not, is colony ii QR but will die next generation?
            } else if (queen_survival[ii] == -1 && population[ii,number_alleles + 4] == 1){
                
                                        #if so, it produces drones but it will not persist in the new population
                
                ## QR colony produces drones
                distribution_alleles[population[ii,1:2]] <- distribution_alleles[population[ii,1:2]] + 
                    (population[ii,number_alleles + 3] * 0.5) 
                
                                        #or is the colony QL?
            } else if (population[ii,number_alleles + 4] == 0){ 
                ## QL colonies produce no new colonies but do produce drones
######
                ## For each QL colony "row", if the QL colony is a "newly-produced daughter colony that became queenless during her mating flight" then the queen nuclear alleles and the spermathecal contents are not from the dead virgin queen but instead from the virgin queen's mother (since we are assuming that the virgin queen is lost on her mating flight, which means that the queenless workers are actually her mother's daughters. In effect, this means that half of the drone alleles from a QL colony come from the original queen (i.e. the virgin queen's mother) since each worker carries one allele from the original queen's nucleus. The other half come from the original queen's spermatheca since each worker carries one allele from the original queen's spermatheca. If the QL colony was headed by a mated queen that died when finding a new nest, then the workers are derived from the mated queen.
                ## The differences between the two QL colonies are accounted for elsewhere (i.e. below all QL colonies are treated identically, since each row contains the relevant info (the virgin queen's mother's nuclear/spermatheca for newly QL colonies but the mated queen's nucleus/spermatheca for established queens that die trying to relocate to a new nest)

                ## original queen's nucleus
                distribution_alleles[population[ii,1:2]] <- distribution_alleles[population[ii,1:2]] + 
                    (population[ii,number_alleles + 3] * 0.25) ## note the 0.25 instead of 0.5 
                ## original queen's spermatheca
                distribution_alleles <- distribution_alleles +  ( (population[ii,number_alleles + 3] * 0.5) * population[ii,3:(number_alleles+2)] ) ## population[ii,3:(number_alleles+2)] is the spermathecal vector that sums to 1, so I need to multiply by this to correctly assign the relative proportion of each. I multiply by 0.5 because the entire thing counts for half (each nuclear allele counts for 0.25, giving 0.5 in total for the nuclear also)
            }
            
                                        #if the iith colony is QR, it can produce daughter swarms (even if it will become dead or QL)
                                        #if a QR colony produces a non-zero number of daughter colonies
            if (population[ii,number_alleles + 4] == 1 && swarming_vector[ii] > 0){
                
                for (jj in 1:swarming_vector[ii]){ #loop through those daughter colonies
                                        #first, choose queen alleles for jjth daughter colony of colony ii
                    
                                        #randomly choose an allele from previous queen (equal chance of either allele)
                    new_population[swarm_counter,1] <- population[ii,sample(1:2,1)]
                    
                                        #randomly choose an allele from the previous queen's spermatheca (probabilty depends on frequency of sperm in spermatheca)
                    new_population[swarm_counter,2] <- sample(1:number_alleles,1, prob = population[ii,3:(number_alleles + 2)])
                    while (new_population[swarm_counter,1] == new_population[swarm_counter,2]){ #if the new queen is homozygous, rechoose
                        new_population[swarm_counter,1] <- population[ii,sample(1:2,1)] 
                        new_population[swarm_counter,2] <- sample(1:number_alleles,1, prob = population[ii,3:(number_alleles + 2)])
                    }
                    ## store which colony this entry is derived from
                    population_original_queen_index[swarm_counter] <- ii
                    
                    swarm_counter <- swarm_counter + 1
                    
                }
            } #if it's QL, it can't produce daughter swarms 
        } ## for ends the for 1:N loop
        
        distribution_alleles <- distribution_alleles / sum(distribution_alleles)
        
                                        #choose which drones each queen mates with 
        sampled_drones <- 
            matrix(sample(1:number_alleles,N_new_generation * number_drone_matings,replace = TRUE, 
                          prob = distribution_alleles), nrow = N_new_generation, ncol = number_drone_matings)
        
                                        #for all NEW colonies ONLY:
        
                                        #transform drone alleles into proportions - drones mate with queens 
                                        #calculate colony fitness
                                        #assign daughter swarms a QR or QL state
        
getColonyFitness <- function(colony_ID, number_alleles, population){

    return( population[colony_ID, number_alleles + 3] )

}

isColonyExtinct <- function(population, number_alleles){
    ## test for extinction (either no colonies added (lhs) or those added are QL (rhs))
    return( !as.logical( NROW(population) ) || !sum( population[ , number_alleles + 4] ) )
}

generateNewPopulation <- function(population, number_alleles, number_drone_matings,
            
            if (sum(new_population[ii,3:(number_alleles + 2)]) == 0){ # means daughter colony NOT surviving queen
                ## now I first chose whether the colony goes QR or QL, as this affects where the workers are derived from

                new_population[ii,number_alleles + 4] <- runif(1) < (1 - probability_QL_colony) # 1 if QR, 0 if QL

                if (new_population[ii,number_alleles + 4] == 1){ ## QR
                    
                    new_population[ii,3:(number_alleles + 2)] <- tabulate(bin = sampled_drones[ii,],nbins = number_alleles) / number_drone_matings
                    
                                        # determine proportion of homozygosity by checking each queen allele and multiplying it by the corresponding drone allele
                                        # each queen allele contributes half the total homozygosity and total homozygosity is multipled by its cost
                    new_population[ii,number_alleles + 3] <- 1 - (((0.5 * new_population[ii,new_population[ii,1] + 2]) + 
                                                                   (0.5 * new_population[ii,new_population[ii,2] + 2])) * cost_homozoygosity)
                } else { ## QL
                    ## population[population_original_queen_index[ii], ] gives the vector corresponding to the original queen 
                    new_population[ii,1:(number_alleles + 3)] <- population[population_original_queen_index[ii], 1:(number_alleles + 3)]
                    ## tag as QL colony
                    new_population[ii,number_alleles + 4] <- 0
                }
                
            } # else, do nothing because surviving queens do not mate again
            
        }
        
        population<- new_population
        return(population)
    }
    
}
