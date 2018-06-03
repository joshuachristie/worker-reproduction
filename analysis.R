num_sims <- 10000
num_files <- 48
num_cols <- 12
counter <- 1
num_gen <- 6

path_to_folder <- sprintf('%s/manuscript_data/', getwd())

summary <- data.frame(matrix(numeric(num_files*num_cols),nrow = num_files,ncol = num_cols))

name_of_sims <- c('QL_no_reproduction','QL_reproduction',
                  'QL_reproduction_multiple_founder','QL_no_reproduction_multiple_founder')

dm <- 25

for (name_sim in name_of_sims){
    for (ch in c(5, 50)){
        for (na in c(15)){
            for (as in c(2,4,6)){
                for (pq in c(30,50)){

                    filename <- sprintf('%s%s_na%d_ch%d_as%d_pq%d_dm%d.RData',path_to_folder,name_sim,na,ch,as,pq,dm)

                    if (file.exists(filename)){
                        
                        
                        load(filename)
                        
                        summary[counter, 1] <- name_sim
                        summary[counter, 2] <- ch/100
                        summary[counter, 3] <- na
                        summary[counter, 4] <- as
                        summary[counter, 5] <- pq/100
                        
                        ## mean alleles at generation num_gen
                        mean_counter <- numeric(num_sims)
                        
                        for (ii in 1:num_sims){
                            
                            sim <- results[[ii,1]]
                            
                            mean_counter[ii] <- sum(sim[num_gen + 1,] > 0) 
                            
                        }
                        
                        summary[counter,6] <- mean(mean_counter)
                        summary[counter,7] <- sd(mean_counter)
                        summary[counter,8] <- median(mean_counter)

                        ## proportions of extinctions
                        summary[counter,9] <- sum(unlist(results[,4]))/(num_sims + sum(unlist(results[,4])))
                        
                        ## work out the mean level of inviable males at generation num_gen (6)
                        ## the only relevant pairing is drone[ii] with queen[ii] (everything else
                        ## gives viable female diploids)
                        mean_proportion_diploid_males <- numeric(num_sims)
                        for (ii in 1:num_sims){

                            queen_dist <- unlist(results[[ii,1]][num_gen+1,])
                            drone_dist <- unlist(results[[ii,2]][num_gen+1,])
                            
                            mean_proportion_diploid_males[ii] <- sum(queen_dist * drone_dist)
                            
                        }

                        summary[counter,10] <- mean(mean_proportion_diploid_males)
                        summary[counter,11] <- sd(mean_proportion_diploid_males)
                        summary[counter,12] <- sprintf('%1.2f--%1.2f',range(mean_proportion_diploid_males)[1],range(mean_proportion_diploid_males)[2])
                        
                        counter <- counter + 1
                        
                    }
                }
            }
        }
    }
}

colnames(summary) <- c("sim type", "cost homozygosity",	"alleles", "swarm", "rate of qlness", "mean", "sd", "median", "proportion of extinctions", "DMP at gen 6", "sd", "range")

output_filename <- sprintf('%s/manuscript_data/data_summary.csv', getwd())

write.csv(summary, output_filename, row.names=FALSE)
