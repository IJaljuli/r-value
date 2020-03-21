# R setup
devtools::install_github( "guido-s/meta"    , force=T )
devtools::install_github( "IJaljuli/metarep", force=T )
library(reshape2)
library(doParallel)
library(doRNG)
library(ggplot2)
library(dplyr)
library(metarep)
library(meta)
library(gridExtra)
library(grid)

## General simulation parameters 

RNG_SEED <- 1219
Rounds <- as.integer(1E5)
SIM_NR_CORES <- 62
R.core <- ceiling(Rounds/SIM_NR_CORES)

###############################################################
########################### Part I  ###########################
###############################################################

## Mixture of fixed true signals, analysed by either fixed  or raandom effects.
# mu1 - the increased effects. 
# mu2 - the decreased effects.
signal.mixture <- function(mu1 = 0 , mu2 = 0 , N1 = 0 , N2 = 0 , 
                           N.studies = 4 , n.i = 25 ,
                           equal.sample.size = T , 
                           sigma.i = 1, R = 10 , U = 2,
                           do.fixed = T , do.random = T ,
                           truncation.t = 1 ){
  # Setting sample sizes matrix. 
  if( !equal.sample.size){
    n.i.1 <- c( 20, 208, 24, 190, 58, 36, 51) + 2
    n.i.1 <- c(  n.i.1  , 15 )
    n.i.2 <- c( 20, 119, 22, 185, 29, 51, 47) + 2
    n.i.2 <- c(  n.i.2  , 16 )
    if(N.studies != length(n.i.1)){
      n.i.1 <- sample( x = n.i.1 , size = N.studies, replace = T)
      n.i.2 <- sample( x = n.i.2 , size = N.studies, replace = T)
    }
  }else{
    if(length(n.i) == 1){ n.i <- rep(n.i , N.studies) }
    n.i.1 <- n.i.2 <- n.i 
  }
  n.i.1.mat <- n.i.1 ; n.i.2.mat <-  n.i.2 
  # n.i.1.mat <- matrix(data = n.i.1 , nrow = 1 , ncol = N.studies)#[rep(1,R),]
  # n.i.2.mat <- matrix(data = n.i.2 , nrow = 1 , ncol = N.studies)#[rep(1,R),]
  
  # Setting true signals of the N.studies
  theta.i <- c(rep(0,N.studies-N1-N2) , rep(mu1,N1) , rep(mu2,N2) )
  
  # Tables for saving the results per round ( columns ) 
  Meta_Analysis.fixed <- # fixed-effects meta
    Meta_Analysis.random <- rep(NA,R)  # random-effects meta
  Replicability_Analysis.fixed <- cbind(  u = U , matrix(NA , nrow = length(U) , ncol = R) )  # replicability analysis of meta using fixed-effects t-test
  Replicability_Analysis.random1 <- matrix(NA , nrow = length(truncation.t)*length(U) , ncol = R)
  # Replicability-analysis by Zaykin will be performed for various values of u and truncaion thresholds t (as inserted). 
  # These values will be indixed in columns 1-2 of the matrices Replicability_Analysis.random1 and Replicability_Analysis.random2: 
  Replicability_Analysis.random1 <- cbind( expand.grid( t = truncation.t  ,u =  U) , 
                                           Replicability_Analysis.random1 ) # registering r(u)-values
  Replicability_Analysis.random2 <- Replicability_Analysis.random1 # Consistency-analysis 
  
  for (r in 1:R){
    # Sample studies means 
    ybar.1 <- rnorm(n = N.studies , sd = sigma.i/sqrt(n.i.1.mat) ) + theta.i
    ybar.2 <- rnorm(n = N.studies , sd = sigma.i/sqrt(n.i.2.mat) )
    # Sample studies means and SEs.
    se.1 <- sqrt( rchisq(n = N.studies , df = n.i.1.mat - 1) / n.i.1.mat ) * sigma.i
    se.2 <- sqrt( rchisq(n = N.studies , df = n.i.2.mat - 1) / n.i.2.mat ) * sigma.i
    
    FmetaModel <- RmetaModel <- NULL
    theres.an.error <- F
    # Fitting RE meta-analysis. If the model fails or p-values == 0, the data will be resampled.
    if ( do.random ){
      RmetaModel <- tryCatch( metacont(n.e = n.i.1.mat , mean.e = ybar.1 , sd.e = se.1, 
                                       n.c = n.i.2.mat , mean.c = ybar.2 , sd.c = se.2,
                                       comb.fixed = F  , comb.random = T , hakn = F,
                                       studlab = paste0('S' , 1:N.studies) , sm = 'MD' ) ,
                              error=function(e) e)
      
      theres.an.error <- inherits( RmetaModel , "error") 
      
      if(!theres.an.error){ # chek if there are less than two valid p-values (i.e. != 0 ). 
        valid.pvs <- sum( RmetaModel$pval >= 10^(-39))
        if ( valid.pvs  < 2 ){
          RmetaModel$pval <- replace( RmetaModel$pval , RmetaModel$pval <  10^(-39) , 10^(-39) )
          valid.pvs <- sum(RmetaModel$pval >= 10^(-39))
          theres.an.error <- ( valid.pvs  < 2 )
        }
      }
      
      # To eliminate cases where all p-values are identically 0, we fix them at 10^-39. If this still doesn't work, resample data. 
      
      while( theres.an.error ){
        
        ybar.1 <- rnorm(n = N.studies , sd = sigma.i/sqrt(n.i.1.mat) ) + theta.i
        ybar.2 <- rnorm(n = N.studies , sd = sigma.i/sqrt(n.i.2.mat) )
        
        se.1 <- sqrt( rchisq(n = N.studies , df = n.i.1.mat - 1) / n.i.1.mat ) * sigma.i
        se.2 <- sqrt( rchisq(n = N.studies , df = n.i.2.mat - 1) / n.i.2.mat ) * sigma.i
        
        RmetaModel <- tryCatch( metacont(n.e = n.i.1.mat , mean.e = ybar.1 , sd.e = se.1, 
                                         n.c = n.i.2.mat , mean.c = ybar.2 , sd.c = se.2,
                                         comb.fixed = do.fixed  , comb.random = T , hakn = F,
                                         studlab = paste0('S' , 1:N.studies) , sm = 'MD' ) ,
                                error=function(e) e)
        theres.an.error <- inherits( RmetaModel , "error")
        if(!theres.an.error){
          valid.pvs <- sum(RmetaModel$pval >= 10^(-39))
          if ( valid.pvs  < 2 ){
            RmetaModel$pval <- replace( RmetaModel$pval , RmetaModel$pval <  10^(-39) , 10^(-39) )
            valid.pvs <- sum(RmetaModel$pval >= 10^(-39))
            theres.an.error <- ( valid.pvs  < 2 )
          }
        }
        
      }
      
      Meta_Analysis.random[r] <- RmetaModel$pval.random  # saving meta-analysis p-value. 
      
      # perform replicability-analysis for the various values of u and truncation thresholds t.  
      for ( j in 1:nrow(Replicability_Analysis.random1) ){
        u.j <- Replicability_Analysis.random1[j,'u'] # u used fo the replicability requirement.
        # Perform one-sided replicability analyses: 
        Replicability_Analysis.random1[j,r+2] <- 
          metarep(x = RmetaModel, u = u.j, 
                  t = Replicability_Analysis.random1[j,'t'], alternative = 'greater',
                  report.u.max = F, common.effect = F)$r.value
        
        Replicability_Analysis.random2[j,r+2] <-
          metarep(x = RmetaModel, u = u.j, 
                  t = Replicability_Analysis.random1[j,'t'], alternative = 'less',
                  report.u.max = F, common.effect = F)$r.value
        
        temp.replicability <- 2*min( c(0.5,
                                       Replicability_Analysis.random1[j,r+2] ,
                                       Replicability_Analysis.random2[j,r+2]) )
        
        temp.inconsistency <- ifelse( u.j > 1 , NA ,  
                                      2*min( c(0.5, max( c(
                                        Replicability_Analysis.random1[j,r+2] ,
                                        Replicability_Analysis.random2[j,r+2])  ) ) ) )
        
        Replicability_Analysis.random2[j,r+2] <- temp.inconsistency
        Replicability_Analysis.random1[j,r+2] <- temp.replicability
        
      }
    }
    
    # Perform fixed-effects meta-analysis and the matching replicability analysis  
    if ( do.fixed ){
      
      FmetaModel <- tryCatch( metacont(n.e = n.i.1.mat , mean.e = ybar.1 , sd.e = se.1,
                                       n.c = n.i.2.mat , mean.c = ybar.2 , sd.c = se.2,
                                       comb.fixed = T , comb.random = F , hakn = F,
                                       studlab = paste0('S' , 1:N.studies) , sm = 'MD' ) ,
                              error=function(e) e)
      
      Meta_Analysis.fixed[r] <- FmetaModel$pval.fixed 
      
      # start here
      for ( j in 1:nrow(Replicability_Analysis.fixed) ){
        u.j <- Replicability_Analysis.fixed[j,'u'] # u used fo the replicability requirement.
        Replicability_Analysis.fixed[j , r + 1 ] <- metarep(x = FmetaModel,u = u.j, t = NULL, alternative = 'two-sided',
                                                            report.u.max = F, common.effect = T)$r.value
      }
      
    }
  }
  
  if ( ! is.null(nrow(Replicability_Analysis.random1)) ){ # Turning the results to vector layouts.
    temp1 <- temp2 <-  us  <- tees <- rounds <- c()
    for ( t_i in 1:nrow(Replicability_Analysis.random1)){
      temp1 <- c(temp1 , as.numeric( Replicability_Analysis.random1[t_i,-c(1,2) ] ) )  
      temp2 <- c(temp2 , as.numeric( Replicability_Analysis.random2[t_i,-c(1,2) ] ) )  
      tees <- c( tees , rep( Replicability_Analysis.random1[t_i, 1 ] , R ) )
      us <- c( us , rep( Replicability_Analysis.random1[t_i, 2 ] ,  R ) )
      rounds <- c( rounds , 1:R )
    }
    Replicability_Analysis.random1 <- temp1
    Replicability_Analysis.random2 <- temp2
    remove(temp2)
  }
  
  if ( ! is.null(nrow(Replicability_Analysis.fixed)) ){ # Turning the results to vector layouts.
    temp1 <-  us.fix <- rounds.fix <- c()
    for ( t_i in 1:nrow(Replicability_Analysis.fixed)){
      temp1 <- c(temp1 , as.numeric( Replicability_Analysis.fixed[t_i,-c(1) ] ) )  
      us.fix <- c( us.fix , rep( Replicability_Analysis.fixed[t_i, 1 ] ,  R ) )
      rounds.fix <- c( rounds.fix , 1:R )
    }
    temp1 ->  Replicability_Analysis.fixed
    remove(temp1)
  }
  
  
  rep.a <- length(Meta_Analysis.fixed) 
  rep.b <- length(Meta_Analysis.random)
  rep.c <- length(Replicability_Analysis.fixed) 
  rep.d <- length(Replicability_Analysis.random1) 
  rep.e <- length(Replicability_Analysis.random2)
  
  # Return:
  data.frame(Model = rep(c('Fixed' , 'Random','Fixed' ,'Random','Random'),
                         c( rep.a , rep.b , rep.c , rep.d , rep.e ) ) , 
             Method = rep(c('Meta-Analysis' , 'Replicability Analysis'), 
                          c(rep.a + rep.b , rep.c + rep.d + rep.e )),
             pvalue = c(Meta_Analysis.fixed, 
                        Meta_Analysis.random,
                        Replicability_Analysis.fixed, 
                        Replicability_Analysis.random1, 
                        Replicability_Analysis.random2),
             truncation.t = c( rep(NA , rep.a + rep.b + rep.c ) , tees , tees ),
             U = c( rep(NA , rep.a + rep.b ) , us.fix ,rep(us , 2 ) ),
             Sigma.i = sigma.i, 
             rounds = c( rep( 1:R , 2 ), rounds.fix, rounds , rounds),
             Tested.alternative =  rep( c( 'Replicability', 'Inconsistency'), c(rep.a + rep.b + rep.c + rep.d , rep.e ) )
  )
  
}

# Setting scenarios for simulatins generatong figures 2(simulation ='2A',FE.meta = F), 3(simulation ='2B') and 10(simulation ='2B',FE.meta = T)
scenarios.table <-
  expand.grid( simulation = '2A', number.of.studies = c(4,8,20),
               Sign.consistency = c('Consistent' ,'Inconsistent') , 
               nonnull.studies.per.side = c(1,2),
               jump.length = c( seq(from = 0.1 , to = 1 , by = .01 ), seq(from = 1.01 , to = 2 , by = 0.05 )),
               Sigma.i = 1 , 
               equal.sample.size = c(F,T),
               FE.meta = F ,stringsAsFactors = F) %>% 
  filter( ! ((Sign.consistency ==  'Inconsistent') & (nonnull.studies.per.side == 1)) )

scenarios.table <- rbind( scenarios.table ,
                          expand.grid( simulation = '2B', number.of.studies = c(4,8,20),
                                       Sign.consistency = 'Consistent' , 
                                       nonnull.studies.per.side = 0:10,
                                       jump.length = c(1,2),
                                       Sigma.i = 1 , 
                                       equal.sample.size = c(T,F), FE.meta = F,stringsAsFactors = F ),
                          expand.grid( simulation = '2B', number.of.studies =  c(4,8,20),
                                       Sign.consistency = 'Consistent' , 
                                       nonnull.studies.per.side = 0:10,
                                       jump.length = c(1,2),
                                       Sigma.i = 1 , 
                                       equal.sample.size = c(T,F), FE.meta = c(F,T),stringsAsFactors = F))

scenarios.table <- scenarios.table %>%
  mutate( Sign.consistency = as.character(Sign.consistency) , 
          Tested.alternative = 'two-sided' )

scenarios.table <- scenarios.table %>% 
  filter( nonnull.studies.per.side <= number.of.studies ) %>%
  filter( (simulation == '2A')& (number.of.studies == 8) & equal.sample.size)
## ggtale: data frame for summing simulation results producing the graphs. 

ggtable <- NULL


scenarios.table <-  anti_join(x = scenarios.table,
                              y = ggtable[,c('jump.length','nonnull.studies.per.side','Sign.consistency') ] )

cl <- makeCluster(SIM_NR_CORES)
registerDoParallel(cl)
clusterExport(cl , c())

for( i in 1:nrow(scenarios.table) ){
  print( paste(i , 'out of ', nrow(scenarios.table) ))
  StartTime = Sys.time()
  # nonnull signals:  
  Mu1 = scenarios.table$jump.length[i] # The size of the increased effect
  Mu2 = ifelse( scenarios.table$Sign.consistency[i] =='Consistent' ,
                NA , -2 * Mu1 ) # The size of the decreased effect
  
  N.nonnull1 = scenarios.table$nonnull.studies.per.side[i] # Number of studies sampled from a population with increased effect
  N.nonnull2 = ifelse( scenarios.table$Sign.consistency[i] =='Consistent' ,
                       0 , ceiling( N.nonnull1 / 2) ) # Number of studies sampled from a population with decreased effect
  
  dof <- scenarios.table$FE.meta[i] # perform fixed-effects meta-analysis
  dor <-  T  # perform random-effects meta-analysis
  
  
  if( as.character(scenarios.table$simulation[i]) == '2B') { # 2B <-> simulations for u=2 and alpha.tilde=0.05
    temp.U <-  2
    temp.t <- 0.05
  }else{
    temp.U <- c( 1, 2, 3 )
    temp.t <- c(1 , .5 , .05 )
  }
  
  # perform paralell computations for each scenario:
  all.runs <-  foreach( j = 1:SIM_NR_CORES, .options.RNG = RNG_SEED ,  .combine = rbind ) %dorng% {
    library(metarep);library(meta); library(dplyr)
    
    single.run <-
      signal.mixture(mu1 = Mu1 , mu2 = Mu2 , 
                     N1 = N.nonnull1 , N2 = N.nonnull2,
                     N.studies =  scenarios.table$number.of.studies[i],
                     equal.sample.size = scenarios.table$equal.sample.size[i],
                     sigma.i = scenarios.table$Sigma.i[i] ,
                     R = R.core,
                     do.fixed = dof , do.random = dor ,
                     truncation.t =  temp.t, U = temp.U )
    single.run %>% mutate(rounds = rounds + (j-1)*R.core)
  }
  
  if(dof){
    all.runs.inconsistency <- NULL
  }else{
    all.runs.inconsistency <- all.runs %>% # Inconsistency: significant right- and left-sided r(u=1). 
      filter( U == 1 , Tested.alternative == 'Inconsistency' ,
              Method == 'Replicability Analysis', Model == 'Random') %>% 
      group_by( Method , Model , truncation.t , Tested.alternative , Sigma.i) %>% 
      summarise( rejection.rate = mean( pvalue < 0.05 )) %>% 
      mutate( U = 0 )
    
  }
  
  all.runs <- all.runs %>% filter(Tested.alternative == 'Replicability' ) %>%
    group_by( Method, Model , truncation.t , U  , Tested.alternative , Sigma.i) %>%
    summarise( rejection.rate = mean( pvalue < 0.05) )
  
  all.runs <- rbind( all.runs , all.runs.inconsistency[, colnames(all.runs) ])
  
  
  
  scenario.results <- data.frame( all.runs ,
                                  simulation = scenarios.table$simulation[i],
                                  number.of.studies = scenarios.table$number.of.studies[i],
                                  Sign.consistency = scenarios.table$Sign.consistency[i] , 
                                  nonnull.studies.per.side = scenarios.table$nonnull.studies.per.side[i],
                                  jump.length = scenarios.table$jump.length[i],
                                  equal.sample.size = scenarios.table$equal.sample.size[i],
                                  FE.meta = dof )

  ggtable  <- rbind( ggtable  , scenario.results )
  EndTime = Sys.time()
  print(EndTime-StartTime)
}
stopCluster(cl)

# prepping the results tables for the graphs: 
ggtable <- NULL
ggtable$Method <- factor(ggtable$Method , 
                         levels = c("Replicability Analysis","Meta-Analysis") )
ggtable$Sign.consistency <- as.character(ggtable$Sign.consistency)
ggtable$Tested.alternative <- as.character(ggtable$Tested.alternative)

# Splitting the results: 
ggtableB <- ggtable %>% filter(simulation == '2B') 
ggtable  <- ggtable%>% filter(simulation == '2A')

# Recoding for facet labels & adding columns for titles:
ggtable <- ggtable %>% 
  mutate( rejection.type = recode(nonnull.studies.per.side , 
                                  '1' = 'Single-deviating study' , '2' = 'Two-deviating studies' ), 
          Sign.consistency = recode( Sign.consistency , 
                                     'Consistent'  = 'Consistent Setting'  , "Inconsistent" = 'Inconsistent Setting'  ),
          U = recode( U , '1' = 'u==1'  , '2' = 'u==2' ,'3' = 'u==3' ,'0' = '0' ))

ggtable <- ggtable %>% mutate( titles.rows =  U  ) %>% 
  mutate(titles.rows = replace(titles.rows , Tested.alternative=='Inconsistency' , 'Inconsistency'))

ggtable$titles.rows <- factor( ggtable$titles.rows ,
                               levels = c( 'u==1' , 'u==2', 'u==3' , 'Inconsistency', 'NA') )

ggtable$titles <- "(~mu*~', 0 , ... , 0 ')"
ggtable <- ggtable %>% 
  mutate( titles = replace(titles,
                           (Sign.consistency=="Consistent Setting")&(rejection.type == "Two-deviating studies"),
                           "(~mu*~','~mu*~' , 0 , ... , 0 ')")) %>% 
  mutate( titles = replace(titles,
                           (Sign.consistency=="Inconsistent Setting")&(rejection.type == "Two-deviating studies"),
                           "(~mu*~','~mu*~','~-2*mu*~' , 0 , ... , 0 ')"))
ggtable$titles <- factor(ggtable$titles , 
                         levels = c("(~mu*~', 0 , ... , 0 ')" , 
                                    "(~mu*~','~mu*~' , 0 , ... , 0 ')",
                                    "(~mu*~','~mu*~','~-2*mu*~' , 0 , ... , 0 ')") )

# More explicit reshaping the results table for the figures
rand.data <-  subset(ggtable , (Model=='Random')&(Tested.alternative != 'inconsistency') )
tem <- rand.data[ rand.data$Method == 'Meta-Analysis', ]
t.vals <- unique(rand.data$truncation.t) ; t.vals <- t.vals[!is.na(t.vals)]
u.vals <- unique(rand.data$U) ; u.vals <- u.vals[!is.na(u.vals)]; u.vals <- u.vals[u.vals>0]
rows.titles.text <- unique(rand.data$titles.rows)[-1]
tem.repli <- expand.grid(row.num = 1:nrow(tem) , truncation.t = t.vals , 
                         U =  u.vals , titles.rows = rows.titles.text)
tem <-  tem[tem.repli$row.num , ] 
tem$truncation.t = tem.repli$truncation.t
tem$U = tem.repli$U 
tem$titles.rows = tem.repli$titles.rows

rand.data <-  rbind( tem[colnames(rand.data)] , 
                     subset(rand.data ,  rand.data$Method != 'Meta-Analysis' ) )
rand.data$truncation.t <- as.factor(rand.data$truncation.t)
rand.data$U <- as.factor(rand.data$U)
rand.data$Method <- factor(rand.data$Method , levels = rev(sort(levels(rand.data$Method))) )
rand.data <- rand.data %>% filter(!is.na(rand.data$titles.rows))

(equal.sampleA <- subset(rand.data , equal.sample.size&(number.of.studies==8)) %>%
    mutate( truncation.t = replace(truncation.t , is.na(truncation.t) , 1 ) ) %>%
    ggplot( .   , aes( x = jump.length , y = rejection.rate  , color = truncation.t))+
    geom_line( aes( linetype = Method ) ) + 
    ylab('Rejection Rate') + xlab(expression(mu)) +
    facet_grid( titles.rows ~ titles , drop=T,scales = 'free',labeller=label_parsed ) +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    #  labs(title = "Constant sample size") + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5) ) +
    scale_color_discrete( labels = c('0.05' , '0.5' , 'No truncation' ),
                          name = 'Truncation at t =' ))


equal.sampleA4 <- subset(rand.data , equal.sample.size&(number.of.studies==4)) %>%
    mutate( truncation.t = replace(truncation.t , is.na(truncation.t) , 1 ) ) %>%
    ggplot( .   , aes( x = jump.length , y = rejection.rate  , color = truncation.t))+
    geom_line( aes( linetype = Method ) ) + 
    ylab('Rejection Rate') + xlab(expression(mu)) +
    facet_grid( titles.rows ~ titles   , drop=T,scales = 'free',labeller=label_parsed ) +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    #  labs(title = "Constant sample size") + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
  theme(strip.text = element_text(size = 7), title = element_text(size = 9))+
  scale_color_discrete( labels = c('0.05' , '0.5' , 'No truncation' ),
                          name = 'Truncation at t =' ) +
    ggtitle( 'Equal Sample Sizes, N = 4 Studies')

equal.sampleA20 <- subset(rand.data , equal.sample.size&(number.of.studies==20)) %>%
    mutate( truncation.t = replace(truncation.t , is.na(truncation.t) , 1 ) ) %>%
    ggplot( .   , aes( x = jump.length , y = rejection.rate  , color = truncation.t))+
    geom_line( aes( linetype = Method ) ) + 
    ylab('Rejection Rate') + xlab(expression(mu)) +
    facet_grid( titles.rows ~ titles , drop=T,scales = 'free',labeller=label_parsed ) +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    #  labs(title = "Constant sample size") + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
  theme(strip.text = element_text(size = 7), title = element_text(size = 9))+
  scale_color_discrete( labels = c('0.05' , '0.5' , 'No truncation' ),
                          name = 'Truncation at t =' ) + 
    ggtitle( 'Equal Sample Size, N = 20 Studies')



nonequal.sampleA <- subset(rand.data , (!equal.sample.size)&(number.of.studies==8)) %>%
    mutate( truncation.t = replace(truncation.t , is.na(truncation.t) , 1 ) ) %>%
    ggplot( .   , aes( x = jump.length , y = rejection.rate  , color = truncation.t))+
    geom_line( aes( linetype = Method ) ) + 
    ylab('Rejection Rate') + xlab(expression(mu)) +
    facet_grid( titles.rows ~ titles , drop=T,scales = 'free',labeller=label_parsed ) +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    #  labs(title = "Constant sample size") + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
  theme(strip.text = element_text(size = 7), title = element_text(size = 9))+
  scale_color_discrete( labels = c('0.05' , '0.5' , 'No truncation' ),
                          name = 'Truncation at t =' ) +
  ggtitle( 'Different Group Sizes, N = 8 Studies') 



nonequal.sampleA4 <- subset(rand.data , (!equal.sample.size)&(number.of.studies==4)) %>%
    mutate( truncation.t = replace(truncation.t , is.na(truncation.t) , 1 ) ) %>%
    ggplot( .   , aes( x = jump.length , y = rejection.rate  , color = truncation.t))+
    geom_line( aes( linetype = Method ) ) + 
    ylab('Rejection Rate') + xlab(expression(mu)) +
    facet_grid( titles.rows ~ titles   , drop=T,scales = 'free',labeller=label_parsed ) +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    #  labs(title = "Constant sample size") + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
  theme(strip.text = element_text(size = 7), title = element_text(size = 9))+
  scale_color_discrete( labels = c('0.05' , '0.5' , 'No truncation' ),
                          name = 'Truncation at t =' ) + 
    ggtitle( 'Different Group Sizes, N = 4 Studies')


nonequal.sampleA20 <- subset(rand.data , (!equal.sample.size)&(number.of.studies==20)) %>%
  mutate( truncation.t = replace(truncation.t , is.na(truncation.t) , 1 ) ) %>%
  ggplot( .   , aes( x = jump.length , y = rejection.rate  , color = truncation.t))+
  geom_line( aes( linetype = Method ) ) + 
  ylab('Rejection Rate') + xlab(expression(mu)) +
  facet_grid( titles.rows ~ titles , drop=T,scales = 'free',labeller=label_parsed ) +
  geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
  #  labs(title = "Constant sample size") + 
  scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
  theme(strip.text = element_text(size = 7), title = element_text(size = 9))+
  scale_color_discrete( labels = c('0.05' , '0.5' , 'No truncation' ),
                        name = 'Truncation at t =' ) + 
  ggtitle( 'Different Group Sizes, N = 20 Studies')


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

fixA.leg <- g_legend(nonequal.sampleA + theme (legend.box = 'horizontal', legend.box.just = 'right') +
                       guides(colour = guide_legend(order = 1),
                                               linetype = guide_legend(order = 2)) )

pdf('sim2_all.pdf' , width = 8 ,height = 11)
equal.sampleA
grid.arrange(nonequal.sampleA + theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL)+ theme(axis.text = element_text(size=7)),
             fixA.leg,
             nonequal.sampleA4+ theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL)+ theme(axis.text = element_text(size=7)),
             equal.sampleA4+ theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL)+ theme(axis.text = element_text(size=7)),
             nonequal.sampleA20+ theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL)+ theme(axis.text = element_text(size=7)),
             equal.sampleA20+ theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL) + theme(axis.text = element_text(size=7)),
             ncol=2, top = textGrob('Extensions of Fig. 2' ,gp=gpar(fontsize=14)))
dev.off()


rand.data <-  subset(ggtableB , (Model=='Random')&(Tested.alternative == 'Replicability') )
rand.data$Method <- factor(rand.data$Method , levels = rev(sort(levels(rand.data$Method))) )
rand.data$jump.length <- as.factor(rand.data$jump.length)
colnames(rand.data)[which(colnames(rand.data)=='jump.length')] = 'Signal'

pdf('sim4_all.pdf' , width = 6,height = 4)
equal.sampleB <- subset(rand.data , equal.sample.size&(number.of.studies==8) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .7,end = .15) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
    theme(strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  )+ 
    ggtitle( 'Equal Sample Sizes \nN = 8 Studies')

equal.sampleB4 <- subset(rand.data , equal.sample.size&(number.of.studies==4) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .7,end = .15) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  ) + 
    ggtitle( 'Equal Sample Sizes \nN = 4 Studies')

equal.sampleB20 <- subset(rand.data , equal.sample.size&(number.of.studies==20) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .7,end = .15) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  ) + 
    ggtitle( 'Equal Sample Sizes\n N = 20 Studies')

nonequal.sampleB20 <- 
    subset(rand.data , (!equal.sample.size)&(number.of.studies==20) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .5,end = 0) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  ) + 
    ggtitle( 'Different Group Sizes \nN = 20 Studies')

(nonequal.sampleB8 <- 
    subset(rand.data , (!equal.sample.size)&(number.of.studies==8) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .5,end = 0) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5) ) )


nonequal.sampleB4 <- 
    subset(rand.data , (!equal.sample.size)&(number.of.studies==4) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .5,end = 0) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  ) + 
    ggtitle( 'Different Group Sizes \nN = 4 Studies')



fix.data <-  rbind( 
  ggtableB %>% filter( Model=='Fixed' , Method == "Replicability Analysis" )  ,
  ggtableB %>% filter( Model=='Random', Method == "Meta-Analysis" ) ) 
fix.data <- fix.data %>%  filter( FE.meta , Tested.alternative == 'Replicability')

fix.data$Method <- factor(fix.data$Method , levels = levels(rand.data$Method))
fix.data$jump.length <- as.factor(fix.data$jump.length)
colnames(fix.data)[which(colnames(fix.data)=='jump.length')] = 'Signal'


equal.sampleB.Fix <- subset(fix.data , equal.sample.size&(number.of.studies==8) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .7,end = .15) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested \n (Common effect)' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  ) + 
  ggtitle( 'Equal Sample Sizes \nN = 8 Studies')

equal.sampleB4.Fix <- subset(fix.data , equal.sample.size&(number.of.studies==4) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .7,end = .15) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested  \n (Common effect)' , 'RE Meta-Analysis') ) + 
    theme(strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  ) + 
    ggtitle( 'Equal Sample Sizes \nN = 4 Studies')

equal.sampleB20.Fix <- subset(fix.data , equal.sample.size&(number.of.studies==20) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .7,end = .15) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested  \n (Common effect)' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  ) + 
    ggtitle( 'Equal Sample Sizes \nN = 20 Studies')

nonequal.sampleB20.Fix <- 
    subset(fix.data , (!equal.sample.size)&(number.of.studies==20) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .5,end = 0) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested  \n (Common effect)' , 'RE Meta-Analysis') ) + 
    theme(strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5), 
          title = element_text(size = 9)  ) + 
    ggtitle( 'Different Group Sizes \nN = 20 Studies')

(nonequal.sampleB8.Fix <- 
    subset(fix.data , (!equal.sample.size)&(number.of.studies==8) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .5,end = 0) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested \n (Common effect)' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5) ))

nonequal.sampleB4.Fix <- 
    subset(fix.data , (!equal.sample.size)&(number.of.studies==4) ) %>%
    ggplot( .   , aes( x = nonnull.studies.per.side , y = rejection.rate ) , color = 'black')+
    geom_line( aes( linetype = Method , color = Signal ) ) +
    geom_point( aes( shape = Signal , color = Signal ) ) +
    scale_colour_grey(start = .5,end = 0) +
    ylab('Rejection Rate') + xlab('Number of studies with non-null signal') +
    geom_abline(slope = 0,intercept = .05 , alpha = 0.5)  + 
    scale_linetype_discrete( labels = c('Suggested  \n (Common effect)' , 'RE Meta-Analysis') ) + 
    theme(legend.position = 'bottom',strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 10) , strip.text.y = element_text(size = 9.5) , 
          title = element_text(size = 9) ) + 
    ggtitle( 'Different Group Sizes \nN = 4 Studies')

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


rand.leg <- g_legend(equal.sampleB)
grid.arrange( equal.sampleB + theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL),
              equal.sampleB4 + theme(legend.position = 'none') + xlab(NULL) + ylab(NULL),
              nonequal.sampleB4 + theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL), 
              rand.leg,
              top = textGrob('Extensions of Fig. 3' ,gp=gpar(fontsize=14)))

fix.leg <- g_legend(equal.sampleB4.Fix)
grid.arrange( equal.sampleB.Fix + theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL),
              equal.sampleB4.Fix + theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL),
              nonequal.sampleB4.Fix + theme(legend.position = 'none') + xlab(NULL)+ ylab(NULL), 
              fix.leg, ncol = 2,
              top = textGrob('Extensions of Fig. 10 (Appendix)' ,gp=gpar(fontsize=14)))



dev.off()




# pdf('sim_2_multi_t_uDec12.pdf' , width = 8, height = 6 )
# 
# equal.sampleA
# nonequal.sampleA
# equal.sampleB
# nonequal.sampleB
# 
# dev.off()



###############################################################
########################### Part II ###########################
###############################################################


RE.assumptions <- function( mu = 0 ,N.studies = 4 , n.i = 25 ,
                            equal.sample.size = T , #rep(25 , N.studies),
                            sigma.i = 1,R = 10 , tau = 1 ,
                            truncation.t = 1 ){
  do.random = T
  
  # Setting sample sizes matrix. 
  if( !equal.sample.size){
    n.i.1 <- c( 20, 208, 24, 190, 58, 36, 51) + 2
    n.i.1 <- c(  n.i.1  , 15 )
    n.i.2 <- c( 20, 119, 22, 185, 29, 51, 47) + 2
    n.i.2 <- c(  n.i.2  , 16 )
    if(N.studies != length(n.i.1)){
      n.i.1 <- sample( x = n.i.1 , size = N.studies, replace = T)
      n.i.2 <- sample( x = n.i.2 , size = N.studies, replace = T)
    }
  }else{
    if(length(n.i) == 1){ n.i <- rep(n.i , N.studies) }
    n.i.1 <- n.i.2 <- n.i 
  }
  # n.i.1.mat <- matrix(data = n.i.1 , nrow = 1 , ncol = N.studies)[rep(1,R),]
  # n.i.2.mat <- matrix(data = n.i.2 , nrow = 1 , ncol = N.studies)[rep(1,R),]
  
  # Tables for saving the results per round ( columns ) 
  Meta_Analysis.random <- I2 <- number.of.positive.signals <-  rep(NA,R)
  Meta_Analysis.fixed <- Replicability_Analysis.fixed <-  Meta_Analysis.random  <-  rep(NA,R)
  u.less.mat <- u.greater.mat <- rvalue.1.max <-
    cbind(  t = truncation.t,  matrix(NA,ncol = R , nrow = length(truncation.t)) )
  
  
  for (r in 1:R){
    theres.an.error <- F
    # Setting true signals of the N.studies
    theta.i <- rnorm(n = N.studies , mean = mu , sd = tau)
    
    # Counting the number of positive true signals ( The alternative is correct according to the signals of the sampled-true effects). 
    number.of.positive.signals[r] <-  sum( theta.i > 0 )
    
    # Sample studies means 
    ybar.1 <- rnorm(n = N.studies , sd = sigma.i/sqrt(n.i.1.mat) ) + theta.i
    ybar.2 <- rnorm(n = N.studies , sd = sigma.i/sqrt(n.i.2.mat) )
    
    # Sample studies means and SEs.
    se.1 <- sqrt( rchisq(n = N.studies , df = n.i.1.mat - 1) / n.i.1.mat ) * sigma.i
    se.2 <- sqrt( rchisq(n = N.studies , df = n.i.2.mat - 1) / n.i.2.mat ) * sigma.i
    
    # Fitting RE meta-analysis. If the model fails or p-values == 0, the data will be resampled.
    RmetaModel <- tryCatch( metacont(n.e = n.i.1.mat , mean.e = ybar.1 , sd.e = se.1, 
                                     n.c = n.i.2.mat , mean.c = ybar.2 , sd.c = se.2,
                                     comb.fixed = F , comb.random = T , hakn = F,
                                     studlab = paste0('S' , 1:N.studies) , sm = 'MD'), error=function(e) e)
    theres.an.error <- inherits( RmetaModel , "error")
    
    if(!theres.an.error){
      valid.pvs <- sum(RmetaModel$pval >= 10^(-39))
      if ( valid.pvs  < 2 ){
        RmetaModel$pval <- replace( RmetaModel$pval , RmetaModel$pval <  10^(-39) , 10^(-39) )
        valid.pvs <- sum(RmetaModel$pval >= 10^(-39))
        theres.an.error <- ( valid.pvs  < 2 )
      }
    }
    
    while( theres.an.error ){
      
      ybar.1 <- rnorm(n = N.studies , sd = sigma.i/sqrt(n.i.1.mat) ) + theta.i
      ybar.2 <- rnorm(n = N.studies , sd = sigma.i/sqrt(n.i.2.mat) )
      
      se.1 <- sqrt( rchisq(n = N.studies , df = n.i.1.mat - 1) / n.i.1.mat ) * sigma.i
      se.2 <- sqrt( rchisq(n = N.studies , df = n.i.2.mat - 1) / n.i.2.mat ) * sigma.i
      
      RmetaModel <- tryCatch( metacont(n.e = n.i.1.mat , mean.e = ybar.1 , sd.e = se.1, 
                                       n.c = n.i.2.mat , mean.c = ybar.2 , sd.c = se.2,
                                       comb.fixed = F , comb.random = T , hakn = F,
                                       studlab = paste0('S' , 1:N.studies) , sm = 'MD'), error=function(e) e)
      theres.an.error <- inherits( RmetaModel , "error")
      
      if(!theres.an.error){
        valid.pvs <- sum(RmetaModel$pval >= 10^(-39))
        if ( valid.pvs  < 2 ){
          RmetaModel$pval <- replace( RmetaModel$pval , RmetaModel$pval <  10^(-39) , 10^(-39) )
          valid.pvs <- sum(RmetaModel$pval >= 10^(-39))
          theres.an.error <- ( valid.pvs  < 2 )
        }
      }
    }
    
    Meta_Analysis.random[r] <- RmetaModel$pval.random #p-value of RE meta-analysis
    I2[r]<- RmetaModel$I2 # I2 as an evaluation of data heterogeniety. 
    
    for ( j in 1:nrow(u.less.mat) ){
      # Check for inconsistency (we declare inconsistency if r^L(1) < 0.05 and r^R(1)<0.05 )
      rvalue.1.less <- metarep(x = RmetaModel, u = 1 , t = u.less.mat[j,'t'], common.effect = F, 
                               alternative = 'less', report.u.max = F )$r.value
      rvalue.1.greater <- metarep(x = RmetaModel, u = 1 , t = u.less.mat[j,'t'], common.effect = F, 
                                  alternative = 'greater', report.u.max = F )$r.value
      rvalue.1.max[j,r+1] <- max(c( rvalue.1.greater  , rvalue.1.less ))
      
      # Find left- and right-bound u_max. later this will be used to test H^(u/n) for all u.  
      u_max_bounds <- find_umax(x = RmetaModel , alternative = 'two-sided',
                                t = u.less.mat[j,'t'] ,common.effect = F )$u_max
      u.less.mat[j,r+1] <- u_max_bounds[2]
      u.greater.mat[j,r+1] <- u_max_bounds[3]
    }
    
  }
  
  if ( ! is.null(nrow(u.less.mat)) ){ # Turning the results to vector layouts.
    temp1<- temp2 <- temp3 <-  tees <- rounds <- c()
    for ( t_i in 1:nrow(u.less.mat)){
      temp1 <- c(temp1 , as.numeric( u.less.mat[t_i,-c(1) ] ) )  
      temp2 <- c(temp2 , as.numeric( u.greater.mat[t_i,-c(1) ] ) )
      temp3 <- c(temp3 , as.numeric( rvalue.1.max[t_i,-c(1) ] ) )
      tees <- c( tees , rep( u.less.mat[t_i, 1 ] , R ) )
      rounds <- c(rounds , 1:R )
    }
    u.less.mat <- temp1
    u.greater.mat <- temp2
    rvalue.1.max <- temp3
    remove(temp1);remove(temp2);remove(temp3)
  }
  
  data.frame(Model = rep(c('Random', 'Random'),
                         c(R,length(truncation.t)*R) ) , 
             Method = rep(c('Meta-Analysis' , 'Replicability Analysis'), 
                          c(R,    length(truncation.t)*R) ),
             pvalue = c(Meta_Analysis.random, rep(NA , length(u.less.mat))),
             rvalue.uIS1 = c(rep(NA, R ), rvalue.1.max),
             u.less = c(rep(NA, R), u.less.mat ),
             u.greater = c(rep(NA, R), u.greater.mat ), 
             truncation.t = c( rep(NA , R ) , tees ),
             I2 = c( I2 ,I2[rounds] ) ,
             rounds = c( rep(1:R , 1 ) , rounds ),
             Sigma.i = mean(sigma.i), 
             Tau = tau, 
             number.of.true.positive = c( rep(NA , R ) ,number.of.positive.signals[rounds] )
  )
}

original.tau <- sqrt(0.5)

scenarios.table <- rbind(
  expand.grid( number.of.studies = c(8,20) ,
               jump.length =  seq(from = 0 , to = 2 , by = 0.04 ),
               Sigma.i = 1 , Tau = c(0.25 , 0.5)   ,
               fig7 = F ,  equal.sample.size = c(F,T) ),
  expand.grid( number.of.studies = 4,
               jump.length =  seq(from = 0 , to = 2 , by = 0.04),
               Sigma.i = 1 ,
               Tau = c( 0.5 , 1 )  ,
               fig7 = F ,  equal.sample.size = c(F,T) ) )

scenarios.table  <- scenarios.table  #%>% filter( !equal.sample.size )
ggtable.all.random <- ggtable.all.random.nnull <- NULL 

cl <- makeCluster(SIM_NR_CORES)
registerDoParallel(cl)
clusterExport(cl , c())

for( i in 1 : nrow(scenarios.table)){
  print( paste(i , 'out of ', nrow(scenarios.table) ))
  Mu = scenarios.table$jump.length[i]
  
  # StartTime = Sys.time()
  
  all.runs <-  foreach( j = 1:SIM_NR_CORES, .options.RNG = RNG_SEED ,  .combine = rbind ) %dorng% {
    library(meta);library(metarep);library(dplyr)
    single.run <-
      RE.assumptions(mu = Mu , tau = scenarios.table$Tau[i],
                     N.studies = scenarios.table$number.of.studies[i] ,
                     sigma.i = scenarios.table$Sigma.i[i] ,
                     R = R.core,
                     truncation.t = 0.05,
                     equal.sample.size =  scenarios.table$equal.sample.size[i])
    #all.runs <-
    single.run %>% mutate(rounds = rounds + (j-1)*R.core, 
                          equal.sample.size = scenarios.table$equal.sample.size[i])
  }
  
  
  all.runs.nnull <- all.runs %>%
    filter(( Model == 'Random' )&(Method=='Replicability Analysis')) %>%
    group_by( Method , Model , truncation.t  , Tau , equal.sample.size, Sigma.i ) %>%
    mutate(nif.inconsistency = mean((number.of.true.positive > 0)&(number.of.true.positive<8 ) ,na.rm = T ),
           nif.1 = mean(number.of.true.positive >= 1 ,na.rm = T ),
           nif.2 = mean(number.of.true.positive >= 2 ,na.rm = T ),
           nif.3 = mean(number.of.true.positive >= 3 ,na.rm = T ),
           nif.4 = mean(number.of.true.positive >= 4 ,na.rm = T ),
           nif.5 = mean(number.of.true.positive >= 5 ,na.rm = T ),
           nif.6 = mean(number.of.true.positive >= 6 ,na.rm = T ),
           nif.7 = mean(number.of.true.positive >= 7 ,na.rm = T ),
           nif.8 = mean(number.of.true.positive >= 8 ,na.rm = T ) )
  
  all.runs.nnull <- all.runs.nnull[,-c(3:6,8,9,12)]
  all.runs.nnull <-  melt(data = all.runs.nnull , 
                          id.vars = c("Model","Method","truncation.t","Sigma.i",'Tau','equal.sample.size') ,
                          value.name = 'null.is.false.fraction', variable.name = 'U_null')
  all.runs.nnull$U_null <- as.character(all.runs.nnull$U_null) %>% 
    strsplit(. , split = 'nif.') %>% 
    lapply(X = .,FUN = function(x) x[2] ) %>% unlist(.)
  
  all.runs.nnull <- all.runs.nnull %>%  
    mutate( U_null = paste0('u~=~',U_null) ) %>%
    mutate( U_null = replace(U_null , U_null=='u~=~inconsistency','Inconsistency' ) )
  
  
  all.runs.inconsistency <- all.runs %>%
    filter( (Method == 'Replicability Analysis')&( Model == 'Random')) %>%
    mutate(null.is.false =  as.character((number.of.true.positive>0)&(number.of.true.positive<8)) ) %>%
    group_by( Method , Model , truncation.t  , Tau , equal.sample.size, null.is.false) %>%
    summarise( rejection.rate05 = mean( rvalue.uIS1 < .05,na.rm = T ), 
               rejection.rate20 = mean( rvalue.uIS1 < .20 ,na.rm = T), 
               I2 = mean(I2,na.rm = T) ) %>%
    mutate( U = 0 , Tested.alternative = 'Inconsistency')
  
  all.runs.inconsistency.temp <- all.runs %>%
    filter( (Method == 'Replicability Analysis')&( Model == 'Random')) %>%
    mutate(null.is.false =  'Either' ) %>%
    group_by( Method , Model , truncation.t  , Tau , equal.sample.size, null.is.false) %>%
    summarise( rejection.rate05 = mean( rvalue.uIS1 < .05,na.rm = T ), 
               rejection.rate20 = mean( rvalue.uIS1 < .20 ,na.rm = T), 
               I2 = mean(I2,na.rm = T) ) %>%
    mutate( U = 0 , Tested.alternative = 'Inconsistency')
  
  all.runs.inconsistency <- rbind(all.runs.inconsistency , 
                                  all.runs.inconsistency.temp[names(all.runs.inconsistency)] )
  remove(all.runs.inconsistency.temp)
  
  temp <- expand.grid(U = 1:8 , x = 1:nrow(all.runs))
  all.runs <- cbind(U = temp$U , all.runs[temp$x , ] )
  remove(temp)
  
  # temp.all.runs <- rbind(
  #   all.runs %>% mutate( Tested.alternative = 'Replicability', 
  #                        null.is.false =  ((8 - number.of.true.positive) >= U),
  #                        rejection.rate20 = (u.less >= U) ) ,
  #   all.runs %>% mutate( Tested.alternative = 'Replicability', 
  #                        null.is.false =  ( number.of.true.positive >= U),
  #                        rejection.rate20 = (u.greater >= U) )) %>%
  #   group_by( Method, Model , truncation.t , U ,Tested.alternative , Tau , equal.sample.size , null.is.false ) %>%
  #   summarise( rejection.rate20 = mean( rejection.rate20 & null.is.false ,na.rm = T) ,
  #              I2 = mean(I2,na.rm = T)  ) %>% 
  #   mutate(rejection.rate05 = rejection.rate20 , null.is.false = as.character(null.is.false ))
  
  temp.all.runs <- 
    all.runs %>% mutate( Tested.alternative = 'Replicability', 
                         null.is.false.L =  ((8 - number.of.true.positive) >= U),
                         rejection.rate20.L = (u.less >= U) ,
                         null.is.false.G =  ( number.of.true.positive >= U),
                         rejection.rate20.G = (u.greater >= U) ) %>%
    mutate( rejection.rate20.L = (rejection.rate20.L & null.is.false.L),
            rejection.rate20.G = (rejection.rate20.G & null.is.false.G),
            null.is.false = (null.is.false.L | null.is.false.G) ) %>%
    group_by( Method, Model , truncation.t , U ,Tested.alternative , Tau , equal.sample.size , null.is.false ) %>%
    summarise( rejection.rate20 = mean( rejection.rate20.L | rejection.rate20.G ,na.rm = T) ,
               I2 = mean(I2,na.rm = T)  ) %>% 
    mutate(rejection.rate05 = rejection.rate20 , null.is.false = as.character(null.is.false ))
  
  
  all.runs.temp.meta <- all.runs %>% filter(Method != 'Replicability Analysis') %>%
    mutate( null.is.false = '' ) %>%
    group_by( Method, Model , truncation.t , U , Tau , equal.sample.size , null.is.false ) %>%
    summarise( rejection.rate20 = mean(pvalue < 0.05,na.rm = T) ,
               I2 = mean(I2,na.rm = T) ) %>% 
    mutate(rejection.rate05 = rejection.rate20 , Tested.alternative = 'Replicability' )
  
  all.runs <- all.runs %>% 
    mutate( Tested.alternative = 'Replicability',
            null.is.false = 'Either') %>%
    group_by( Method, Model , truncation.t , U ,Tested.alternative , Tau , equal.sample.size , null.is.false ) %>%
    summarise( rejection.rate20 = mean((u.greater >=U) | (u.less >=U),na.rm = T) ,
               I2 = mean(I2,na.rm = T) ) %>% 
    mutate(rejection.rate05 = rejection.rate20 )
  
  
  all.runs <- rbind( all.runs ,  all.runs.temp.meta[colnames(all.runs)], 
                     temp.all.runs[colnames(all.runs)]  ,all.runs.inconsistency[names(all.runs) ])
  
  remove(temp.all.runs); remove(all.runs.temp.meta); remove(all.runs.inconsistency)
  
  
  
  scenario.results <- data.frame( all.runs ,
                                  number.of.studies = scenarios.table$number.of.studies[i],
                                  jump.length = scenarios.table$jump.length[i],
                                  fig7 =  scenarios.table$fig7[i] )
  
  ggtable.all.random <- rbind( ggtable.all.random, scenario.results )
  
  
  all.runs.nnull <- data.frame( all.runs.nnull ,
                                number.of.studies = scenarios.table$number.of.studies[i],
                                jump.length = scenarios.table$jump.length[i],
                                fig7 =  scenarios.table$fig7[i] )
  ggtable.all.random.nnull <- rbind(ggtable.all.random.nnull , all.runs.nnull )
  
  # EndTime = Sys.time()
  print( scenarios.table[i,] )
  print(EndTime - StartTime)
  
}
stopCluster(cl)

ggtable.all.random <- NULL

ggtable.all.random <- ggtable.all.random %>% filter(null.is.false != 'FALSE')
ggtable.all.random$null.is.false <- replace(ggtable.all.random$null.is.false,
                                            ggtable.all.random$null.is.false == 'TRUE' , 'Null~is~false')
ggtable.all.random$null.is.false <- replace(ggtable.all.random$null.is.false,
                                            ggtable.all.random$null.is.false == 'Either' , 'Null~is~either~ture~or~false')

ggplot(ggtable.all.random, aes(x=I2)) +
  geom_histogram() + facet_grid( number.of.studies+  equal.sample.size ~ Tau ) +
  geom_vline(xintercept = c(0.5,0.7))

ggtable.all.random <- ggtable.all.random %>%
  mutate(I2.size = 'Moderate~heterogeneity' , sample.size = 'Equal~sample~size') %>% 
  mutate( I2.size = replace( I2.size , (Tau == 1)|( (Tau==0.5)&(number.of.studies>4) ) , 'High~heterogeneity'),
          sample.size = replace(sample.size , !equal.sample.size , 'Unequal~sample~size' ))

ggtable.all.random <- ggtable.all.random %>% filter( Model == 'Random' ) #%>% mutate(fig7 = F)
ggtable.all.random$Method <- factor(ggtable.all.random$Method ,
                                    levels = c( 'Replicability Analysis','Meta-Analysis') )
ggtable.all.random$truncation.t <- as.factor(ggtable.all.random$truncation.t )
ggtable.all.random$Tau <- round( ggtable.all.random$Tau , digits = 2 )

# Duplicate the meta-analysis results for all truncarion values. 
rand.data <-  subset(ggtable.all.random,  Method != 'Meta-Analysis' )
tem <- subset(ggtable.all.random,  Method == 'Meta-Analysis' )
if(nrow(tem) == 0 ){ tem <- matrix(tem , nrow = 1) }
t.vals <- unique(rand.data$truncation.t) ; t.vals <- t.vals[!is.na(t.vals)]
u.vals <- unique(rand.data$U) ; u.vals <- u.vals[!is.na(u.vals)] ;  u.vals <- u.vals[u.vals>0]
null.false <- c('Null~is~false','Null~is~either~ture~or~false')
tem.repli <- expand.grid(row.num = 1:nrow(tem) , t =  t.vals , U = u.vals , null.is.false = null.false)
tem <-  tem[tem.repli$row.num , ] 
tem$truncation.t = tem.repli$t
tem$U = tem.repli$U 
tem$null.is.false = tem.repli$null.is.false 
rand.data <-  rbind( tem[colnames(rand.data)] ,  rand.data )

rand.data <- rand.data %>% mutate( U.text=paste0( "u*~'= ", U,"'" ))
rand.data$Tau <-  paste0( "~tau*~'= ", round(rand.data$Tau,digits = 2),"'" )
rand.data$Tested.alternative <- factor(rand.data$Tested.alternative , levels = rev(unique(rand.data$Tested.alternative)) )

ttt <- rand.data #subset( rand.data , (!fig7)&( desired.I2 == .5 ) )
ttt$U.text <- replace(ttt$U.text , ttt$Method == 'Meta-Analysis' , "u*~'= 1'" )
ttt$U.text <- replace(ttt$U.text , ttt$Tested.alternative == 'Inconsistency' , NA )
U.labels <- paste( paste( 'expression(' , 
                          c(unique( ttt$U.text )[-c(9)], 'Inconsistency' ) , ')' ) , collapse = ',' )
U.labels <- paste0('c(' , U.labels ,')')
# ttt$U.text <- replace(ttt$U.text , ttt$U == 0 , NA )
ttt$Tested.alternative <- factor(ttt$Tested.alternative ,
                                 levels = rev(levels(ttt$Tested.alternative)))

ttt <- ttt %>% mutate( desired.I2.text = paste0(   "I^2 == ~", 100*I2 ),
                       equal.sample.size = as.character(equal.sample.size))%>%
  filter(!is.na(rejection.rate20))

ttt$null.is.false <- factor(ttt$null.is.false , levels = c( 'Null~is~either~ture~or~false' , 'Null~is~false' ) )


sim3A_equal <- 
  ttt %>% filter(null.is.false == 'Null~is~either~ture~or~false', number.of.studies == 8, !as.logical(equal.sample.size)) %>% 
  ggplot(  . , 
           aes( x = jump.length , y = rejection.rate05, color = U.text ))+
  geom_line(aes( linetype = Method) ,alpha = 1, lwd = .8) +
  ylab('Rejection Rate') + xlab(expression(mu)) + 
  geom_abline(slope = 0 , intercept =  .05 , color = 'grey')+
  scale_y_discrete(limits=c(0,0.05, 0.5, 1 ),
                   expand = expand_scale(add = c(.08, 0.1 )),
                   labels=c('0',"5%", '50%', "100%" ))+
  facet_grid(   . ~ I2.size     ,labeller = label_parsed , drop = T ) + 
  scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') )+
  scale_color_discrete(labels = eval(parse(text = U.labels ) ) ,
                       name = element_blank() )+
  theme(legend.position = 'right',strip.text = element_text(size = 9),title = element_text(size = 8),
        legend.text.align = 0) + ggtitle('Equal Sample Sizes, N=8 Studies')



sim3A_nonequal <- 
  ttt %>% filter(null.is.false == 'Null~is~either~ture~or~false', number.of.studies == 8, !as.logical(equal.sample.size)) %>% 
  ggplot(  . , 
           aes( x = jump.length , y = rejection.rate05, color = U.text ))+
  geom_line(aes( linetype = Method) ,alpha = 1, lwd = .8) +
  ylab('Rejection Rate') + xlab(expression(mu)) + 
  geom_abline(slope = 0 , intercept =  .05 , color = 'grey')+
  scale_y_discrete(limits=c(0,0.05, 0.5, 1 ),
                   expand = expand_scale(add = c(.08, 0.1 )),
                   labels=c('0',"5%", '50%', "100%" ))+
  facet_grid(   . ~ I2.size     ,labeller = label_parsed , drop = T ) + 
  scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') )+
  scale_color_discrete(labels = eval(parse(text = U.labels ) ) ,
                       name = element_blank() )+
  theme(legend.position = 'right',strip.text = element_text(size = 11),
        legend.text.align = 0) #+ ggtitle('Different Group sizes')
# dev.off()
# 
# pdf('sim3B_all.pdf')

sim3B_equal <- 
  ttt %>% filter(null.is.false == 'Null~is~either~ture~or~false' , number.of.studies != 8,!as.logical(equal.sample.size)) %>% 
  mutate(number.of.studies = factor(paste0(number.of.studies, '~Studies'),
                                     levels=c('4~Studies','20~Studies')) ) %>%
  ggplot(  . , 
           aes( x = jump.length , y = rejection.rate05, color = U.text ))+
  geom_line(aes( linetype = Method) ,alpha = 1, lwd = .8) +
  ylab('Rejection Rate') + xlab(expression(mu)) + 
  geom_abline(slope = 0 , intercept =  .05 , color = 'grey')+
  scale_y_discrete(limits=c(0,0.05, 0.5, 1 ),
                   expand = expand_scale(add = c(.08, 0.1 )),
                   labels=c('0',"5%", '50%', "100%" ))+
  facet_grid(   number.of.studies ~ I2.size     ,labeller = label_parsed , drop = T ) + 
  scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') )+
  scale_color_discrete(labels = eval(parse(text = U.labels ) ) ,
                       name = element_blank() )+
  theme(legend.position = 'right',strip.text = element_text(size = 9),title = element_text(size = 8),
        legend.text.align = 0) +  ggtitle('Equal Sample Sizes')

sim3B_nonequal <- 
  ttt %>% filter(null.is.false == 'Null~is~either~ture~or~false', number.of.studies != 8, !as.logical(equal.sample.size)) %>% 
  mutate(number.of.studies = factor(paste0(number.of.studies, '~Studies'),
                   levels=c('4~Studies','20~Studies')) ) %>%
  ggplot(  . , 
           aes( x = jump.length , y = rejection.rate05, color = U.text ))+
  geom_line(aes( linetype = Method) ,alpha = 1, lwd = .8) +
  ylab('Rejection Rate') + xlab(expression(mu)) + 
  geom_abline(slope = 0 , intercept =  .05 , color = 'grey')+
  scale_y_discrete(limits=c(0,0.05, 0.5, 1 ),
                   expand = expand_scale(add = c(.08, 0.1 )),
                   labels=c('0',"5%", '50%', "100%" ))+
  facet_grid(   number.of.studies ~ I2.size     ,labeller = label_parsed , drop = T ) + 
  scale_linetype_discrete( labels = c('Suggested' , 'RE Meta-Analysis') )+
  scale_color_discrete(labels = eval(parse(text = U.labels ) ) ,
                       name = element_blank() )+
  theme(legend.position = 'right',strip.text = element_text(size = 9), 
        legend.box = 'horizontal', legend.box.just = 'right',
        title = element_text(size = 8),
        legend.text.align = 0) + ggtitle('Different Group sizes')

sim3_leg <- g_legend(sim3B_nonequal + theme(legend.text = element_text(size=7),
                                            legend.key.size = unit(.5, "cm")) +
                       guides(colour = guide_legend(order = 2),
                              linetype = guide_legend(order = 1)))

pdf('sim3A_all.pdf', height = 5, width = 9)
sim3A_nonequal 
dev.off()

pdf('sim3B_all.pdf', height = 7, width = 10)

grid.arrange(sim3A_equal + theme(legend.position = 'none',plot.margin=unit(c(t=.2,r = .8,b = 0, l = 0),"cm")),
             sim3_leg,
             sim3B_nonequal+ theme(legend.position = 'none'),
             sim3B_equal+ theme(legend.position = 'none'),
             ncol = 2, heights=c(1,1.5),
             top = textGrob('Extensions of Fig. 4' ,gp=gpar(fontsize=14)))


dev.off()















