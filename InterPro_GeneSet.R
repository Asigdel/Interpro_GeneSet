## ===========================================================================
## ipDataSource
# 2019-09-20
## ---------------------------------------------------------------------------
## A class to store mapping information between genes, pathways and
## interPro domains
## ---------------------------------------------------------------------------

setClass("ipDataSource",
         representation(genes="character", pathways="character",
                        domains="character", gene2Domains="list",
                        path2Domains="list", dims="numeric",
                        type="character"))

## ===========================================================================
# Step 2
## ===========================================================================

# For all pathways in an 'ipDataSource' object, compute chisquare-based
## similarity of a gene list to each pathway and provide subsampling
## distributions and p-values. Output of the function is a
## list with items:
##    dist: a list containing the subsampling distributions
##    similarity: a list of distance measures of the gene
##                  list to each pathway
##    pvalue: a named vector of p-values indicating similarity of
##            gene list to each pathway
compSimilarities <- function(geneList, dataSource, n=10000, verbose=TRUE){
  geneList <- geneList[geneList %in% dataSource@genes]
  np <- length(dataSource@pathways)
  dists <- setDists <- list()
  ## setup progress report
  if(verbose){
    require(prada)
    mess <- paste("resampling gene list\nwith ", n ,"samples")
    sub <- paste("(pathway 1 of ", np, ")", sep="")
    progress(message=mess, sub=sub)
    on.exit(killProgress())
  }
  ## all domains of the geneList
  listDoms <- unique(unlist(dataSource@gene2Domains[geneList], 
                            use.names=FALSE))
  ## iterate over all pathways
  for(i in seq(along=dataSource@pathways)){
    p <- dataSource@pathways[i]
    dists[[p]] <- resampleGeneLists(pathway=p,
                                    len=length(geneList), n=n,
                                    dataSource=dataSource)
    setDists[[p]] <- sim2Pathway(pathway=p, doms=listDoms,
                                 dataSource=dataSource)
    ## report progress
    if(verbose)
      updateProgress((i)/np*100, autoKill=TRUE,
                     sub= paste("(", i, " of ", np, ")\n",p, sep=""))
  }
  ## compute pvalues from sample distributions
  pvals <- mapply(function(x,y) sum(x>y)/n, dists, setDists)
  return(list(dist=dists, similarity=setDists, pvalue=pvals))
  
  ## ===========================================================================
  # Step 3
  ## ===========================================================================
  
  
  
  dataSource <- function(mapping, type="generic"){
    if(is.null(names(mapping)))
      stop("'mapping' must be named list of pathway mappings to genes\n")
    pathways <- as.character(names(mapping))
    genes <- as.character(unique(unlist(mapping)))
    
    ## get corresponding interpro domains from biomaRt
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    tmp <- getBM(attributes=c("entrezgene", "interpro"), filters="entrezgene",
                 values=genes, mart = ensembl)
    gene2Domains <- split(tmp$interpro, tmp$entrezgene, drop=FALSE)
    missing <- setdiff(genes, names(gene2Domains))
    gene2Domains[missing] <- ""
    domains <- unique(unlist(gene2Domains))
    domains <- domains[!is.na(domains)]
    path2Domains <- lapply(mapping, function(x, gene2Domains)
      unique(unlist(gene2Domains[x])), gene2Domains)
    
    ## the lengths of pathway, gene and domain vectors
    dims <- c(pathway=length(pathways), gene=length(genes),
              domain=length(domains))
    
    
    ## create an object of class 'ipDataSource'
    return(new("ipDataSource", genes=genes, pathways=pathways,
               domains=domains, gene2Domains=gene2Domains,
               path2Domains=path2Domains, dims=dims, type=type))
  }
  
  ## ===========================================================================
  # Step 4
  ## ===========================================================================
  
  ## a wrapper function that fetches the available KEGG data and
  ## corresponding InterPro IDs from biomaRt for the universe of
  ## entrezgene IDs
  
  getKEGGdata <- function(universe){
    ## get all available human KEGG pathway annotations and
    ## corresponding entrezgene IDs
    hKEGGids <- grep("^hsa", ls(KEGGPATHID2EXTID), value=TRUE)
    path2Genes <- mget(hKEGGids, KEGGPATHID2EXTID)
    hKEGGgenes <- union(universe, unique(unlist(path2Genes, use.names=FALSE)))
    hKEGGgenes <-  hKEGGgenes[!is.na(hKEGGgenes)]
    
    ## get corresponding interpro domains from biomaRt
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    tmp <- getBM(attributes=c("entrezgene", "interpro"), filters="entrezgene",
                 values=hKEGGgenes, mart = ensembl)
    gene2Domains <- split(tmp$interpro, tmp$entrezgene, drop=FALSE)
    missing <- setdiff(hKEGGgenes, names(gene2Domains))
    gene2Domains[missing] <- ""
    hKEGGdomains <- unique(unlist(gene2Domains))
    hKEGGdomains <- hKEGGdomains[!is.na(hKEGGdomains)]
    path2Domains <- lapply(path2Genes, function(x, gene2Domains)
      unique(unlist(gene2Domains[x], use.names=FALSE)), 
      gene2Domains)
    
    ## the lengths of pathway, gene and domain vectors
    dims <- c(pathway=length(hKEGGids), gene=length(hKEGGgenes),
              domain=length(hKEGGdomains))
    
    ## create an object of class 'ipDataSource'
    return(new("ipDataSource", genes=hKEGGgenes, pathways=hKEGGids,
               domains=hKEGGdomains, gene2Domains=gene2Domains,
               path2Domains=path2Domains, dims=dims, type="KEGG"))
  }
  
  ## ===========================================================================
  # Step 5
  ## ===========================================================================
  
  ## Map KEGG Ids to description of pathways. Essentially, this is a wrapper
  ## around the KEGGPATHID2NAME environment.
  
  getKEGGdescription <- function(ids)
    data.frame(ID=ids, description=unlist(mget(substring(ids, 4,100),
                                               KEGGPATHID2NAME)))
  ## ===========================================================================
  # Step 6
  ## ===========================================================================
  
  ## The main function to run geneset enrichment based on Interpro
  ## domain signatures. It takes an ipDataSource object and the gene list
  ## as entrezgene IDs.
  
  gseDomain <- function(dataSource, geneset, n=10000, verbose=TRUE){
    ## validate arguments
    missing <- setdiff(geneset, dataSource@genes)
    if(length(missing)>0)
      stop("There are genes missing from the universe:\n",
           paste("   ", missing,collapse="\n"))
    
    ## start computation
    res <- compSimilarities(geneset, dataSource, n=n, verbose=verbose)
    return(res)
  }
  
  ## ===========================================================================
  # Step 7
  ## ===========================================================================
  ## ==========================================================================
  ## The show method for 'ipDataSource' objects
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  setMethod("show", signature="ipDataSource",
            definition=function(object){
              cat("object of class 'ipDataSource' containing mapping data",
                  "\nfor type '", object@type, "' pathways\n", sep="")
              cat("   ", length(object@pathways), " pathways with " ,
                  length(object@genes), " genes\n   and ", length(object@domains),
                  " annotated InterPro domains\n", sep="")
            }
  )
  
  ## ===========================================================================
  # Step 8
  ## =========================================================================== 
  ## Resample n random geneLists of length 'len' from all
  ## universe genes and compute chisquare-based distance measures
  ## to signature of 'pathway'
  
  resampleGeneLists <- function(pathway, len, n, dataSource){
    ## the unique pathway domains
    #ps <- get(pathway, dataSource@path2Domains)
    ps <- dataSource@path2Domains[[pathway]]
    lps <- length(ps)
    ## number of unique domains in universe (cumsum of cont table)
    lu <- length(dataSource@domains)
    sGenes <- replicate(n, sample(seq(along=dataSource@genes), len), 
                        simplify=FALSE)
    gs <- lapply(sGenes, function(x,y) unique(unlist(y[x], use.names=FALSE)),
                 dataSource@gene2Domains)
    x <- sapply(gs, function(x,y) sum(x%in% y), ps)
    y <- lps-x
    n1 <- listLen(gs)
    n2 <- lu-n1
    pval <- sage.test(x,y,n1,n2)
    return(-log10(pval))
  }
  ## ===========================================================================
  # Step 9
  ## ===========================================================================  
  
  
  ## A vectorized version of the modified chisquare statistic in package
  ## sagenhaft. 'x' and 'y' are parallel entries in the contingency table,
  ## 'n1' and 'n2' are the respective outer sums:
  ##      |   |   |
  ##    --|-------|--
  ##      | x | y | 
  ##      |...|...|
  ##    --|---|---|--
  ##      |n1 | n2|
  
  sage.test <- function (x, y, n1 = sum(x), n2 = sum(y)){
    if (any(is.na(x)) || any(is.na(y)))
      stop("missing values not allowed")
    x <- round(x)
    y <- round(y)
    if (any(x < 0) || any(y < 0))
      stop("x and y must be non-negative")
    if (length(x) != length(y))
      stop("x and y must have same length")
    n1 <- round(n1)
    n2 <- round(n2)
    if(length(n1)==1) n1 <- rep(n1, length(x))
    if(length(n2)==1) n2 <- rep(n2, length(x))
    if (length(n1) != length(n2) | length(n1) != length(x))
      stop("n1 and n2 must have same length as x and y")
    if (!missing(n1) && any(x > n1))
      stop("x cannot be greater than n1")
    if (!missing(n2) && any(y > n2))
      stop("y cannot be greater than n2")
    size <- x + y
    p.value <- rep(1, length(x))
    if (any(n1 == n2)) {
      i <- size > 0 & n1 == n2
      if (any(i)) {
        xI <- pmin(x[i], y[i])
        sizeI <- size[i]
        p.value[i] <- pbinom(xI, size = sizeI, prob = 0.5) +
          pbinom(sizeI - xI + 0.5, size = sizeI, prob = 0.5,
                 lower.tail = FALSE)
      }
    }
    if (any(n1 != n2)) {
      prob <- n1/(n1 + n2)
      if (any(size > 10000 & n1 != n2)) {
        big <- size > 10000 & n1 != n2
        ibig <- (1:length(x))[big]
        for (i in ibig)
          p.value[i] <- chisq.test(matrix(c(x[i], y[i], n1[i] - x[i], n2[i] - y[i]), 2, 2))$p.value
      }
      size0 <- size[size > 0 & size <= 10000 & n1 != n2]
      prob0 <- prob[size > 0 & size <= 10000 & n1 != n2]
      mar0 <- unique(cbind(size0, prob0), MARGIN=1)
      if (nrow(mar0))
        for (ind in 1:nrow(mar0)) {
          isize <- mar0[ind,1]
          iprob <- mar0[ind,2]
          i <- (size == isize) & (prob==iprob) & n1 != n2
          p <- dbinom(0:isize, p = (prob[i])[1], size = isize)
          o <- order(p)
          cumsump <- cumsum(p[o])[order(o)]
          p.value[i] <- cumsump[x[i] + 1]
        }
    }
    p.value[p.value>1] <- 1
    return(p.value)
  }
  
  ## ===========================================================================
  # Step 10
  ## ===========================================================================  
  
  
  
  ## Compute the chisquare-based similarity of a gene list to
  ## a pathway
  
  sim2Pathway <- function(pathway, doms, dataSource){
    ## the unique pathway domains
    ps <- dataSource@path2Domains[[pathway]]
    lps <- length(ps)
    ## number of unique domains in universe (cumsum of cont table)
    lu <- dataSource@dims[3]
    ## the contingency table
    x <- sum(doms %in% ps)
    y <- lps-x
    n1 <- length(doms)
    n2 <- lu - n1
    ## modified chisquare test  
    pval <- sage.test(x,y,n1,n2)
    return(-log10(pval))
  }
  
  ## ===========================================================================
  # Step 11
  ## ===========================================================================  
  
  .onLoad <- function(lib, pkg){
    require(methods)
  }
  
  
  