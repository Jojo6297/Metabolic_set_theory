
# This script supplements the manuscript "Leveraging agent-based modeling and co-occurrence data to
# validate a generalized metabolic model of species interaction in the human gut microbiome" by 
# Jyoti Jyoti, Hannah Zoller, Wolfgang zu Castell, and Marc-Thorsten Hütt.

# This script takes as input
# - the genome-scale models of all 73 species (sorted by random set in the folders "random_set_models_[]")
# - the minimal media of all 73 species (sorted by random set in the folders "random_set_mm_[]")

# This script produces the Rdata files "simulations_random_set_[]_repetition_[].Rdata"

# Written by Hannah Zoller, 2025


#Set working directory

work_folder <- "work_folder"
setwd(work_folder)

#Load packages

library(xlsx)
library(stringr)
library(BacArena)
library(doParallel)
library(parallel)
library(matrixStats)
library(ggplot2)
library(Hmisc)
library(parallel)

#Set the number of the random set (1-10) and the number of the repetition (1-10)

random_set_nr <- 
repetition_nr <- 

#Import the genome-scale models 
  
models_vmh <- list()
files <- list.files(paste("/random_set_models_",random_set_nr,sep=""), full.names = TRUE)

for(i in 1:25){
  models_vmh[[i]] <- readMATmod(files[[i]])
}

species_names <- rep(NA,25)

names <- list.files(paste("/random_set_models_",random_set_nr,sep=""))
for(i in 1:length(names)){
  name <- names[[i]]
  species_names[i] <- paste(str_sub(name, start = 1, end = -5), ", AGORA version 1.03",sep="")
}

#Import the minimal media 

mm <- list()

files <- list.files(paste("/random_set_mm_",random_set_nr,sep=""), full.names = TRUE)

for(i in 1:25){
  mm[[i]] <- (read.table(files[[i]], sep=",", header = TRUE))[,1]
}

#Determine the complete medium

arena <- Arena(n = 20,m = 20)
for(s in 1:25){
  model <- models_vmh[[s]]
  bac <- Bac(model=model)
  arena <- addOrg(arena, bac, amount=1)
}
complete_medium <- arena@mediac
c <- length(complete_medium)

richness <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

#Set up the cluster

cl <- makeCluster(10, type="PSOCK",outfile=paste("simulations_random_set_",random_set_nr,"_repetition_",repetition_nr,".log",sep=""))
clusterEvalQ(cl, library(BacArena))
clusterEvalQ(cl,library(matrixStats))
clusterExport(cl, "models_vmh")
clusterExport(cl, "mm")
clusterExport(cl, "complete_medium")
clusterExport(cl, "c")
clusterExport(cl, "richness")
clusterExport(cl, "species_names")
clusterExport(cl, "feeding_amount")

print(system.time(sim_list <- parLapply(cl, 1:10, function(rn){

  list_tupel <- list()
  
  r <- 0
  
  for(a in 1:24){
    for(b in (a+1):25){
      if(a !=b){
        
        sub <- sample(1:c, floor(richness[rn]*c), replace = FALSE)
        medium_frac <- complete_medium[sub]
        medium <- unique(c(mm[[a]],mm[[b]],medium_frac))
        
        r <- r+1
        
        list_singles <- list()
        
        arena <- Arena(n = 10,m = 10)
        model <- models_vmh[[a]]
        bac <- Bac(model=model)
        arena <- addOrg(arena, bac, amount=5)
        model <- models_vmh[[b]]
        bac <- Bac(model=model)
        arena <- addOrg(arena, bac, amount=5)
        
        arena <- addSubs(arena, mediac = medium, smax = 0.1, unit = "fmol/cell")
        sim <- simEnv(arena, time=10)
        abund <- plotGrowthCurve(sim, ret_data = TRUE,use_biomass = TRUE)
        abund_penult <- abund[(abund$time)==1,3]
        abund_ult <- abund[(abund$time)==10+1,3]
        growthrate_comb <- abund_ult/abund_penult
        
        abund_vec <- rep(NA,10)
        abund <- plotGrowthCurve(sim, ret_data = TRUE)
        for(ab in 1:10){
          abund_vec[ab] <- sum(abund[((ab-1)*2+1):((ab-1)*2+2),3])
        }
        
        #Assemble lists of cross-feeding
        
        spec <- c(a,b)
        cons <- list()
        cons[[1]] <- c(1)
        cons[[2]] <- c(1)
        flux <- list()
        flux[[1]] <- c(1)
        flux[[2]] <- c(1)
        time <- list()
        time[[1]] <- c(1)
        time[[2]] <- c(1)
        
        et_alive <- max(which(abund_vec!=0))
        
        for(step in 1:et_alive){
          ff <- findFeeding3(sim, step, mets = complete_medium, plot = FALSE)
          if(length(ff)>0){
            if(length(ff)>0){
              for(sp in 1:2){
                species <- species_names[spec[sp]]
                rows <- which(ff[,2]==species)
                if(length(rows)>0){
                  cons_vec <- cons[[sp]]
                  cons_vec <- c(cons_vec,as.matrix(ff[rows,3]))
                  cons[[sp]] <- cons_vec
                  flux_vec <- flux[[sp]]
                  flux_vec <- c(flux_vec,rowMins(abs(as.matrix(ff[rows,5:6]))))
                  flux[[sp]] <- flux_vec
                  time_vec <- time[[sp]]
                  time_vec <- c(time_vec,rep(step,length(rows)))
                  time[[sp]] <- time_vec
                }
              }
            }
          }
        }
        
        for(step in (et_alive+1):10){
          for(sp in 1:2){
            species <- species_names[spec[sp]]
            cons_vec <- cons[[sp]]
            cons_vec <- c(cons_vec,NA)
            cons[[sp]] <- cons_vec
            flux_vec <- flux[[sp]]
            flux_vec <- c(flux_vec,NA)
            flux[[sp]] <- flux_vec
            time_vec <- time[[sp]]
            time_vec <- c(time_vec,NA)
            time[[sp]] <- time_vec
          }
        }

        for(sp in 1:2){
          cons_vec <- cons[[sp]]
          cons_vec <- cons_vec[-1]
          cons[[sp]] <- cons_vec
          flux_vec <- flux[[sp]]
          flux_vec <- flux_vec[-1]
          flux[[sp]] <- flux_vec
          time_vec <- time[[sp]]
          time_vec <- time_vec[-1]
          time[[sp]] <- time_vec
        }
        
        list_singles[[1]] <- cons
        list_singles[[2]] <- flux
        list_singles[[3]] <- time
        
        #Assemble list of total influx
        
        cons_total <- list()
        cons_total[[1]] <- c(1)
        cons_total[[2]] <- c(1)
        flux_total <- list()
        flux_total[[1]] <- c(1)
        flux_total[[2]] <- c(1)
        time_total <- list()
        time_total[[1]] <- c(1)
        time_total[[2]] <- c(1)
        
        pp <- plotSpecActivity(sim, subs = complete_medium, ret_data = TRUE)
        if(dim(pp)[1]>0){
          for(sp in 1:2){
            species <- species_names[spec[sp]]
            for(step in 1:et_alive){
              r1 <- which(pp$time==step)
              r2 <- which(pp$spec==species)
              r3 <- which(pp$mflux<0)
              rows_intersect1 <- intersect(r1,r2)
              rows_intersect <- intersect(rows_intersect1, r3)
              if(length(rows_intersect)>0){
                cons_total_vec <- cons_total[[sp]]
                cons_total_vec <- c(cons_total_vec,as.matrix(pp[rows_intersect,2]))
                cons_total[[sp]] <- cons_total_vec
                flux_total_vec <- flux_total[[sp]]
                flux_total_vec <- c(flux_total_vec,abs(as.matrix(pp[rows_intersect,3])))
                flux_total[[sp]] <- flux_total_vec
                time_total_vec <- time_total[[sp]]
                time_total_vec <- c(time_total_vec,rep(step,length(rows_intersect)))
                time_total[[sp]] <- time_total_vec
              }
            }
          }
        }
        
        for(sp in 1:2){
          species <- species_names[spec[sp]]
          for(step in 1:et_alive){
            cons_total_vec <- cons_total[[sp]]
            cons_total_vec <- c(cons_total_vec,NA)
            cons_total[[sp]] <- cons_total_vec
            flux_total_vec <- flux_total[[sp]]
            flux_total_vec <- c(flux_total_vec,NA)
            flux_total[[sp]] <- flux_total_vec
            time_total_vec <- time_total[[sp]]
            time_total_vec <- c(time_total_vec,NA)
            time_total[[sp]] <- time_total_vec
          }
        }

        for(sp in 1:2){
          cons_total_vec <- cons_total[[sp]]
          cons_total_vec <- cons_total_vec[-1]
          cons_total[[sp]] <- cons_total_vec
          flux_total_vec <- flux_total[[sp]]
          flux_total_vec <- flux_total_vec[-1]
          flux_total[[sp]] <- flux_total_vec
          time_total_vec <- time_total[[sp]]
          time_total_vec <- time_total_vec[-1]
          time_total[[sp]] <- time_total_vec
        }
        
        list_singles[[4]] <- cons_total
        list_singles[[5]] <- flux_total
        list_singles[[6]] <- time_total
        
        list_singles[[7]] <- medium
        list_singles[[8]] <- abund_vec
        
        rm(sim)
        
        list_tupel[[r]] <- list_singles
        
      }
    }
  }
  
  return(list_tupel)
  
}) ))

save(sim_list, file = paste("simulations_random_set_",random_set_nr,"_repetition_",repetition_nr,".Rdata", sep = ""))
