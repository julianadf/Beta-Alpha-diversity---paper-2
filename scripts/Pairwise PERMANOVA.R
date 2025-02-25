
#'@author Pedro Martinez Arbizu 
#'
#'@examples
#' data(iris)
#' pairwise.adonis2(iris[,1:4]~Species,data=iris)
#'
#' #For strata (blocks), Jari Oksanen recommends in the help of adonis2 to define the 
#' #permutation matrix outside the adonis2 call
#' #In this example I have adapted the adonis2 example to have 3 factors in NO3
#'
#' dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10,30)),field=gl(3,1) )
#' Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(18)/2
#' Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(18)/2
#' Y <- data.frame(Agropyron, Schizachyrium)
#' perm <- how(nperm = 199)
#' setBlocks(perm) <- with(dat, field)
#' adonis2(Y ~ NO3, data = dat, permutations = perm)
#'
#'
#' # pairwise comparison
#' pairwise.adonis2(Y ~ NO3, data = dat, strata = 'field') #notice the apostrophes in strata = 'field' !
#'
#' #this will give same results a doing adonis2 pairwise one by one
#'
#' #for factors '0' and '10'
#' dat2 <- dat[dat$NO3 %in% c('0','10'),]
#' Y2 <- Y[dat$NO3 %in% c('0','10'),]
#' setBlocks(perm) <- with(dat2, field)
#' adonis2(Y2 ~ NO3, data = dat2, permutations = perm)
#' # and so on...
#'
#'@export pairwise.adonis2
#'@importFrom utils combn
#'@importFrom vegan adonis vegdist


pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- ad[1]
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 


