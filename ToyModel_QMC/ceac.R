ceac <- function(N,m){
    pkg <- c("Matrix","MASS","boot","tictoc","sna","ggplot2","RColorBrewer","foreach","reshape2")
    lapply(pkg, library, character.only = TRUE)
    
    tic()
    set.seed(666)
    
    mean_rec <- c(0.99,1.33)
    var_rec <- matrix(c(0.22,0.15,0.15,0.2),2,2)
    mean_rel <- c(-1.48,-0.4)
    var_rel <- matrix(c(0.14,0.05,0.05,0.11),2,2)
    
    lor_rec <- mvrnorm(N,mean_rec,var_rec)
    lor_rel <- mvrnorm(N,mean_rel,var_rel)
    P_rec <- matrix(0,N,3)
    P_rel <- matrix(0,N,3)
    P_rec[,1] <- rbeta(N,6,200)
    P_rec[,2] <- inv.logit(logit(P_rec[,1])+lor_rec[,1])
    P_rec[,3] <- inv.logit(logit(P_rec[,1])+lor_rec[,2])
    
    P_rel[,1] <- rbeta(N,2,100)
    P_rel[,2] <- inv.logit(logit(P_rel[,1])+lor_rel[,1])
    P_rel[,3] <- inv.logit(logit(P_rel[,1])+lor_rel[,2])
    
    C_rec <- rnorm(N,1000,50)
    C_rel <- rnorm(N,2000,100)
    C_norec <- rnorm(N,2500,125)
    
    Q_rec <- rnorm(N,26,2)
    Q_rel  <- rnorm(N,23,3)
    Q_norec <- rnorm(N,20,4)
    
    C_t <- t(replicate(N,c(0,300,30)))
    lambda <- seq(1e3, 5e4, length.out = m)
    CEAC <- matrix(0,m,3)
    
    for(i in seq(1, m)){
        QALY <- P_rec*(1-P_rel)*Q_rec + P_rec*P_rel*Q_rel + (1-P_rec)*Q_norec
        Cost <- P_rec*(1-P_rel)*C_rec + P_rec*P_rel*C_rel + (1-P_rec)*C_norec + C_t
        
        NB <- lambda[i]*QALY - Cost
        NBmax <- replicate(3,apply(NB,1,max))
        CEAC[i,] <- colSums(NB == NBmax)/N
    }
    
    toc()
    
    dt <- data.frame("lambda" = lambda, "NoTreatment" = CEAC[,1], "CBT" = CEAC[,2], "Antidepressant" = CEAC[,3])
    dt <- melt(dt, id="lambda") 
    names(dt) <- c("lambda",'methods','yValue')
   
    my_palette <- c("NoTreatment" = "tomato", "CBT" = "black", "Antidepressant" = "#2166ac")
                    
    p <- ggplot(data = dt) + scale_color_manual(values = my_palette) + 
        geom_line(aes(x = lambda, y = yValue, color = methods)) + ylab("CEAC")+ xlab(expression(lambda)) 
    p
}