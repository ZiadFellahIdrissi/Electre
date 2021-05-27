Electre_1 <-
  function(performanceMatrix,alternatives,criteria,criteriaWeights, minmaxcriteria,
           concordance_threshold=1,discordance_threshold=0){
    
    cat("\014")
    
    if (is.null(dim(performanceMatrix))) 
      stop("less than 2 criteria or 2 alternatives")    
    
    ## check the input data
    
    if (!((is.matrix(performanceMatrix) || (is.data.frame(performanceMatrix))))) 
      stop("wrong performanceMatrix, should be a matrix or a data frame")
    
    if (!(is.vector(alternatives)))
      stop("alternatives should be a vector")
    if (!is.character(alternatives))	
      stop("alternatives should be a character vector")
    if(!(nrow(performanceMatrix)==length(alternatives)))
      stop("length of alternatives should be checked") 	
    
    
    if (!(is.vector(criteria)))
      stop("criteria should be a vector")
    if (!is.character(criteria))	
      stop("criteria should be a character vector") 
    if(!(ncol(performanceMatrix)==length(criteria)))
      stop("length of criteria should be checked") 	
    
    if (!(is.vector(minmaxcriteria)))
      stop("minmaxcriteria should be a vector")
    minmaxcriteria=tolower(minmaxcriteria)
    if(!(ncol(performanceMatrix)==length(minmaxcriteria)))
      stop("length of minmaxcriteria should be checked") 
    n=length(minmaxcriteria)	
    
    for (i in 1:n){
      if(!((minmaxcriteria[i]=='min') ||(minmaxcriteria[i]=='max'))){
        stop(" Vector minmaxcriteria must contain 'max' or 'min' ")
      }
    }
    
    if (!(is.vector(criteriaWeights)))
      stop("criteriaWeights should be a vector")
    if (!is.numeric(criteriaWeights))	
      stop("criteriaWeights should be a numeric vector") 
    if(!(ncol(performanceMatrix)==length(criteriaWeights)))
      stop("length of criteriaWeights should be checked") 
    
    ##     *************************************  End of checking the validity of the "inputs" **********************************************************####    
    
    
    #  *************************************************************   Variables transformation *********************************************************   #	 
    pm  	<- performanceMatrix
    av  	<-  as.list(alternatives)
    cv  	<- as.list(criteria)
    vp  	<- criteriaWeights
    mmv  <- minmaxcriteria
    s1 	<- concordance_threshold
    s2 	<- discordance_threshold
    
    if (!(is.na(match("min", mmv)))) {
      for(j in 1:ncol(pm)){
        valmax=max(pm[,j])
        if (mmv[j]=="min") {
          for (i in 1:nrow(pm)){
            pm[i,j]=valmax-pm[i,j]
          }
        }
      }
      
    }
    
    
    #  ****************************************           Calcul matrice des concordance                ***********************************   #	  
    c1=0
    c2=0
    
    i=1
    n=nrow(pm)
    mc <- matrix (rep(0, n*n), n, n)
    diag(mc)=1
    for (i in 1:nrow(pm)){
      
      k=i+1
      while (k<=n){
        for (j in 1:ncol(pm)){
          
          tp=pm[i,j]-pm[k,j]
          if (tp>= 0) {
            c1=c1+vp[j]
          }
          if (tp<= 0) {
            c2=c2+vp[j]
          }
        }
        mc[k,i]=c1/sum(vp)
        mc[i,k]=c2/sum(vp)
        c1=0
        c2=0
        k=k+1
      }
      
      concordance=mc
    }
    print("La matrice des indices de concordance")
    print(concordance)
  
    
    #  *******************************                Calcul matrice des discordances              ***************************************   #
    
    max1=0
    max2=0
    
    md <- matrix (rep(0, n*n), n, n)
    
    sigma=0						# sigma is  the maximum difference  gi(b)-gi(a) on any criterion i
    for (s in 1:ncol(pm)){
      sigtemp=max(pm[,s])-min(pm[,s])
      if (sigtemp>= sigma){
        sigma=sigtemp
      }
    }
    
    for (i in 1:nrow(pm)){
      
      k=i+1
      while (k<=n){
        for (j in 1:ncol(pm)){	  
          tp1=pm[k,j]-pm[i,j]
          if (tp1>= max1) {
            max1=tp1
          }
        }
        for (j in 1:ncol(pm)){	
          tp2=pm[i,j]-pm[k,j]
          if (tp2>= max2) {
            max2=tp2
          }
        }
        md[k,i]=max1/sigma
        md[i,k]=max2/sigma
        max1=0
        max2=0
        k=k+1
      }
      
      discordance =md
    }
    cat("\n")
    print("La matrice des indices de discordance")
    print(discordance)
    
    
    matfiltre <- matrix (rep(0, n*n), n, n)
    
    for (i in 1:n){
      for (j in 1:n){
        if (i != j){
          if((mc[i,j]>=s1)&&(md[i,j]<=s2)){
            matfiltre[i,j]=1
          }
        }
      }
      
      filtre=matfiltre
    }
    
    
    rownames(matfiltre)=av
    colnames(matfiltre)=av
    
    matgraph=t(matfiltre)
    g1<-graph.adjacency(matgraph); 
    g2=plot(g1)
    
    
    mc <- round(matrix(mc,ncol=n,nrow=n),digits=4)
    md <- round(matrix(md,ncol=n,nrow=n),digits=4)    
    rownames(pm)=av
    colnames(pm)=cv
    rownames(mc)=av
    colnames(mc)=av
    rownames(md)=av
    colnames(md)=av
    rownames(matfiltre)=av
    colnames(matfiltre)=av
    
    
    
    # prepare the output
    sink("result.txt")
    cat("------------------------------------------------------------------------------  Performance table    -------------------------------------------------------------------------","\n")
    cat(" ","\n")
    cat(" ","\n")
    print(pm)
    cat(" ","\n")
    cat(" ","\n")
    cat("--------------------------------------------------------------------------  Concordance  matrix    ---------------------------------------------------------------------------","\n")
    cat(" ","\n")
    cat(" ","\n")
    print(t(mc))
    cat(" ","\n")
    cat(" ","\n")
    cat("---------------------------------------------------------------------------  Discordance  matrix    ---------------------------------------------------------------------------","\n")
    cat(" ","\n")
    cat(" ","\n")
    print(t(md))
    cat(" ","\n")
    cat(" ","\n")
    cat("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------","\n")
    sink()
    
    mc=t(mc)
    md=t(md)  
    out <-  list( "Performance Matrix"=pm,"Concordance Matrix"=mc,"Discordance Matrix" =md)       
  }



get_performanceMatrix <- function(ma){
  
  n=nrow(ma)
  c1 = c()
  c2 = c()
  c3 = c()
  pm <- matrix(data = 0, nrow = n, ncol = 3)
  for (i in 1:n){
    c1 <- append(c1,mean(ma[,i]))
    c2 <- append(c2,min(ma[,i][ma[,i]>7]))
    c3 <- append(c3,var(ma[,i]))
  }
  pm[,1] = c1
  pm[,2] = c2
  pm[,3] = c3
  
  return(pm) 
  
}
help(matrix)


notes_etduiants <-cbind(
  c(7, 13, 8, 12, 11),
  c(8, 11, 11, 12, 11),
  c(20, 2, 10, 3, 15),
  c(16, 14, 16, 14, 13),
  c(12, 12, 8, 8, 10)
)



performanceMatrix <- get_performanceMatrix(notes_etduiants)




alternatives <-c("e1","e2","e3","e4","e5")
criteria <-c("cr1","cr2","cr3")


criteriaWeights <- c(0.5,0.25,0.25)

minmaxcriteria <-c("max","max","min")

Electre_1(performanceMatrix,
          alternatives,
          criteria,
          criteriaWeights,
          minmaxcriteria,
          concordance_threshold=0.4,discordance_threshold=0.3)





