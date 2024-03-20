# betalinkINIS function

# Author: Carmelo Gomez-Martinez
# Last update: 21/03/2024
# contact: carmelogzmz@gmail.com
# github: https://github.com/carmelogzmz/betalinkINIS

# This function "betalinkINIS" is a modification of betalinkr function from
# Frund (Frund 2021; included in bipartite package, Dormann 2021), to include
# independent data on community composition to be considered in the estimation
# of the components of interaction turnover.

# betalinkINIS function calculates interaction turnover and its components
# (rewiring and species turnover) form Poisot's and Novotny's methods, and their
# qualitative and quantitative versions, just in the same way betalink function
# (Poisot) and betalinkr function (Frund) do. However, betalinkINIS allows the
# inclusion of independent data on species community composition (one or both
# trophic levels), that together with the interaction sampling allow to
# calculate the rewiring and species turnover components considering the
# species that occur in the communities although they might not be interacting
# in one of the compared networks.

# betalinkINIS function accepts the following parameters:
# Exclusive of betalinkINIS:
# - n1 and n2: temporal interaction networks to be compared (objects of 
#              class "matrix").
# - c1 and c2: list of species in the communities at time 1 and 2, respectively.
#              each of them (c1 or c2) must be an object of class "list" 
#              containing two vectors: one vector with the names of the species 
#              of the lower trophic level and another with the names of the 
#              species of the higher trophic level.
#              (i.e. c1 = list(c("Plant1", "Plant2", "Plant3"), c("Insect1, "Insect2")))
# Shared with betalinkr function (an extended description of the following
# components can be found in the vignette of the bipartite package):
# - index: the dissimilarity index used. For Poisot's method, it can be one from
#          functions "vegdist" or "betadiver" (from vegan package). For 
#          Novotny's method, has to be 'sorensen' (= 'bray') or 'jaccard'.
# - partitioning: 'poisot' to calculate the Poisot's components of interaction
#                 turnover, or 'cdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==ommondenom' (as is called in betalinkr, Frund 
#                 2021) to calculate the Novotny's components.
# - binary: to indicate if using qualitative or quantitative data.
# - proportions: to standardize the data to proportions prior calculate the
#                components. It is automatically set to FALSE if binary = TRUE.
# - function.dist: function to calculate dissimilarity: 'vegdist' or 'betadiver'
#                  from vegan package.
# - distofempty: if there are no links due to rewiring or there are only 
#                rewiring links, 'zero' return 0 and 'na' return NA.
# - partition.st: partitioning of the species turnover component to discern
#                 between he dissimilarity due to the absence of  species of the
#                 lower trophic level and the one due to the absence of species 
#                 of the higher trophic level.
#                 Only applicable for 'partitioning="commondenom"'
# - partition.rr: partitioning into "true" turnover components (WN.repl, 
#                 OS.repl), and a component due to richness diferences (WN.rich)
#                 Only applicable to 'partitioning="commondenom"'

# CODE
betalinkINIS <- function(n1, n2, c1 = NULL, c2 = NULL, index = "bray", binary = TRUE, partitioning = "commondenom", proportions = !binary, function.dist = "vegdist", distofempty = "zero", partition.st = FALSE, partition.rr = FALSE) {

  require(betalink)
  require(tidyverse)
  require(vegan)

  # 01. QUALITATIVE OR QUANTITATIVE MATRICES
  ## it standardizes to proportions for calculating quantitative estimates
  if (proportions) {
    if (binary) {
      warning("standardizing to proportions for binary index; do you really want to do this?!?")
    }
    n1.sums <- n1
    n2.sums <- n2
    for(row in 1:nrow(n1.sums)) {
      n1.sums[row,] <- rep(sum(n1), ncol(n1.sums))
    }
    for(row in 1:nrow(n2.sums)) {
      n2.sums[row,] <- rep(sum(n2), ncol(n2.sums))
    }
    n1 <- n1 / n1.sums
    n2 <- n2 / n2.sums
  }

  ## Binarize matrix if binary = TRUE
  if(binary) {
    n1[n1 > 0] <- 1
    n2[n2 > 0] <- 1
  }

  # 02. LISTING SPECIES IN NETWORKS
  ## extract species list (v1 and v2) from each network (n1 and n2) and list shared network species (vs)
  v1 <- c(rownames(n1), colnames(n1))
  v2 <- c(rownames(n2), colnames(n2))
  vs <- v1[v1 %in% v2]

  # build communities (ct1 and ct2) from network data (n1, n2) plus independent data (c1, c2), and list shared community species (cts)
  if(!is.null(c1) && !is.null(c2)) {
    ct1.l <- list(sort(c(c1[[1]][!c1[[1]] %in% rownames(n1)], rownames(n1))), sort(c(c1[[2]][!c1[[2]] %in% colnames(n1)], colnames(n1))))
    ct2.l <- list(sort(c(c2[[1]][!c2[[1]] %in% rownames(n2)], rownames(n2))), sort(c(c2[[2]][!c2[[2]] %in% colnames(n2)], colnames(n2))))
    ct1 <- c(ct1.l[[1]], ct1.l[[2]])
    ct2 <- c(ct2.l[[1]], ct2.l[[2]])
    cts <- sort(ct1[ct1 %in% ct2])
  }

  # 03. BUILDING MATRICES FOR TIME 1 AND 2 INCLUDING IN EACH ALL THE SPECIES IN n1 AND n2. 
  pla <- sort(unique(c(rownames(n1), rownames(n2))))
  ins <- sort(unique(c(colnames(n1), colnames(n2))))

  ## create matrix with all the species in n1 and n2
  nt <- matrix(0, nrow = length(pla), ncol = length(ins), dimnames = list(pla, ins))

  ## nt1 and nt2 contain the same interactions than n1 and n2, respectively, but have each all the species in n1 and n2
  filling <- function(n, ncomm) {
    for(row in 1:nrow(ncomm)) {
      for(col in 1:ncol(ncomm)) {
        if(rownames(ncomm)[row] %in% rownames(n) && colnames(ncomm)[col] %in% colnames(n)) {
          prow <- which(rownames(n) == rownames(ncomm)[row])
          pcol <- which(colnames(n) == colnames(ncomm)[col])
          ncomm[row,col] <- n[prow, pcol]
        }
      }
    }
    return(ncomm)
  }

  nt1 <- filling(n1, nt)
  nt2 <- filling(n2, nt)

  ## remove species with no interactions in any of the networks (it should not be any, unless you work with simulated networks)
  cleaning <- function(ncomm1, ncomm2) {
    for(row in nrow(ncomm1):1) {
      if((sum(ncomm1[row,]) + sum(ncomm2[row,])) == 0) {
        ncomm1 <- ncomm1[-row,]
        ncomm2 <- ncomm2[-row,]
      }
    }
    for(col in ncol(ncomm1):1) {
      if((sum(ncomm1[,col]) + sum(ncomm2[,col])) == 0) {
        ncomm1 <- ncomm1[,-col]
        ncomm2 <- ncomm2[,-col]
      }
    }
    return(list(ncomm1, ncomm2))
  }

  nt.list <- cleaning(nt1, nt2)

  nt1 <- nt.list[[1]]
  nt2 <- nt.list[[2]]

  # "matrix2links" reorganize the 2 networks (nt1 and nt2) in 1, so that plant-pollinator species pairs are in columns and networks ID (nt1 and nt2) in rows
  matrix2links <- function(web1, web2) {
    l <- NULL
    v1 <- NULL
    v2 <- NULL
    for(row in 1:nrow(web1)) {
      for(col in 1:ncol(web1)) {
        l <- c(l, paste0(rownames(web1)[row],"__",colnames(web1)[col]))
        v1 <- c(v1, web1[row,col])
        v2 <- c(v2, web2[row,col])
      }
    }
    linkmxm <- matrix(rbind(v1,v2), nrow = 2, ncol = length(l))
    rownames(linkmxm) <- c("comm1","comm2")
    colnames(linkmxm) <- l
    return(linkmxm)
  }
  
  linknt <- matrix2links(nt1, nt2) # complete networks
  
  # 04. BUILDING MATRICES ONLY WITH INTERACTIONS OF SPECIES SHARED BETWEEN NETWORKS
  ## the matrices conserved the dimensions of nt1 and nt2, but only the links involving species shared between networks are filled. 
  nt1.sh.net <- nt1
  nt2.sh.net <- nt2

  # if both matrices have interactions in the same row, keep the values of the cells, otherwise fill the row with 0.
  for(row in 1:nrow(nt1.sh.net)) {
    if((length(nt1[row,][nt1[row,] > 0]) > 0) != (length(nt2[row,][nt2[row,] > 0]) > 0)) {
      nt1.sh.net[row,] <- rep(0, ncol(nt1.sh.net))
      nt2.sh.net[row,] <- rep(0, ncol(nt2.sh.net))
    }
  }

  # if both matrices have interactions in the same column, keep the values of the cells, otherwise fill the col with 0.
  for(col in 1:ncol(nt1.sh.net)) {
    if((length(nt1[,col][nt1[,col] > 0]) > 0) != (length(nt2[,col][nt2[,col] > 0]) > 0)) {
      nt1.sh.net[,col] <- rep(0, nrow(nt1.sh.net))
      nt2.sh.net[,col] <- rep(0, nrow(nt2.sh.net))
    }
  }

  # Using again "matrix2links" function to merge the two networks in one.
  # fill it with the data con nt1.sh.net and nt2.sh.net
  linknt.sharedsp.net <- matrix2links(nt1.sh.net, nt2.sh.net) # links of shared interacting species
        
  # 05. BUILDING MATRICES ONLY WITH INTERACTIONS OF SPECIES SHARED BETWEEN COMMUNITIES
  ## the matrices conserve the dimensions fo nt1 and nt2, but only the links involving species shared between communities are filled. 
  if(!is.null(c1) || !is.null(c2)) {
    nt1.sh.com <- nt1
    nt2.sh.com <- nt2

    # if both matrices have interactions in the same row, keep the values of the cells, otherwise fill the row with 0.
    for(row in 1:nrow(nt1.sh.com)) {
      if((length(nt1[row,][nt1[row,] > 0]) > 0) != (length(nt2[row,][nt2[row,] > 0]) > 0)) {
        if(!rownames(nt1)[row] %in% cts || !rownames(nt2)[row] %in% cts) {
          nt1.sh.com[row,] <- rep(0, ncol(nt1.sh.com))
          nt2.sh.com[row,] <- rep(0, ncol(nt2.sh.com))
        }
      }
    }

    # if both matrices have interactions in the same column, keep the values of the cells, otherwise fill the col with 0.
    for(col in 1:ncol(nt1.sh.com)) {
      if((length(nt1[,col][nt1[,col] > 0]) > 0) != (length(nt2[,col][nt2[,col] > 0]) > 0)) {
        if(!colnames(nt1)[col] %in% cts || !colnames(nt2)[col] %in% cts) {
          nt1.sh.com[,col] <- rep(0, nrow(nt1.sh.com))
          nt2.sh.com[,col] <- rep(0, nrow(nt2.sh.com))
        }
      }
    }
    # Using again "matrix2links" function to merge the two networks in one.
    # fill it with the data con nt1.sh.com and nt2.sh.com
    linknt.sharedsp.com <- matrix2links(nt1.sh.com, nt2.sh.com) # links of community shared species
  }

  # 06. BUILDING MATRICES ADDING THE SPECIES IN THE COMMUNITIES TO THE NETWORKS
  # ncomm1 will have the species and interactions of nt1 plus the species that interact in nt2 and are in ct1 (list of species in the community at time 1). 
  # ncomm2 is build similarly.
  
  if(!is.null(c1) && !is.null(c2)) {
    plat <- sort(unique(c(ct1.l[[1]], ct2.l[[1]])))
    inst <- sort(unique(c(ct1.l[[2]], ct2.l[[2]])))

    ncomm <- matrix(0, nrow = length(plat), ncol = length(inst))
    rownames(ncomm) <- plat
    colnames(ncomm) <- inst

    ncomm1 <- ncomm
    ncomm2 <- ncomm

    for(row in 1:nrow(ncomm)) {
      for(col in 1:ncol(ncomm)) {
        if(rownames(ncomm1)[row] %in% rownames(nt1) && colnames(ncomm1)[col] %in% colnames(nt1)) {
          prow <- which(rownames(nt1) == rownames(ncomm1)[row])
          pcol <- which(colnames(nt1) == colnames(ncomm1)[col])
          ncomm1[row, col] <- nt1[prow, pcol]
        }
        if(rownames(ncomm2)[row] %in% rownames(nt2) && colnames(ncomm2)[col] %in% colnames(nt2)) {
          prow <- which(rownames(nt2) == rownames(ncomm2)[row])
          pcol <- which(colnames(nt2) == colnames(ncomm2)[col])
          ncomm2[row, col] <- nt2[prow, pcol]
        }
      }
    }
  }

  specnt.all <- function(web1, web2) {
    specnt.l1 <- rowSums(web1)
    specnt.l2 <- rowSums(web2)
    specnt.h1 <- colSums(web1)
    specnt.h2 <- colSums(web2)

    specnt.lower <- rbind(specnt.l1, specnt.l2)
    rownames(specnt.lower) <- c("comm1","comm2")

    specnt.higher <- rbind(specnt.h1, specnt.h2)
    rownames(specnt.higher) <- c("comm1","comm2")

    specnt.higher.unique <- specnt.higher[, !(colnames(specnt.higher) %in% colnames(specnt.lower)), drop = F]
    if (is.null(colnames(specnt.higher))) {
      specnt.higher.unique <- specnt.higher
    }
    specnt.all <- cbind(specnt.lower, specnt.higher.unique)
    duplicolnames <- setdiff(colnames(specnt.higher), colnames(specnt.higher.unique))
    specnt.all[, duplicolnames] <- specnt.all[, duplicolnames] + specnt.higher[, duplicolnames]

    return(specnt.all)
  }

  specnt.all.net <- specnt.all(nt1, nt2)

  if(!is.null(c1) && !is.null(c2)) {
    specnt.all.com <- specnt.all(ncomm1, ncomm2)

    for(col in 1:ncol(specnt.all.com)) {
      if(specnt.all.com[1, col] == 0 && colnames(specnt.all.com)[col] %in% ct1) {
        specnt.all.com[1, col] <- 0.1
      }
      if(specnt.all.com[2, col] == 0 && colnames(specnt.all.com)[col] %in% ct2) {
        specnt.all.com[2, col] <- 0.1
      }
    }
  }

  # 07. RUNNING DIFFERENT METHODS
  
  # POISOT METHOD (Poisot et al. 2012)
  if(partitioning == "poisot") {
    if(partition.st || partition.rr) {
      warning("further partitioning only available with method partitioning = 'commondenom'")
    }

    poisotMethod <- function(linknt.sharedsp, specnt.all) {
      # conserved links
      linknt.sharedli <- linknt
      linknt.sharedli[, colSums(linknt.sharedli > 0) == 1] <- 0
      linknt.rewiring <- linknt.sharedsp - linknt.sharedli

      # changes due to rewiring
      linknt.RewSha <- linknt.sharedsp
      # changes due to species turnover
      linknt.uniquesp <- linknt - linknt.sharedsp
      # links conserved + changes due to turnover
      linknt.UniSha <- linknt.uniquesp + linknt.sharedli

      if (proportions) {
        linknt.RewSha <- decostand(linknt.RewSha, method = "total")
        linknt.UniSha <- decostand(linknt.UniSha, method = "total")
      }

      # calculate dissimilarity of species, (b_s), interactions (b_wn) and rewiring (b_os) (by vegdist)
      if (function.dist == "vegdist") {
        b_s <- vegdist(specnt.all, method = index, binary = binary)
        b_wn <- vegdist(linknt, method = index, binary = binary)
        b_zero <- b_wn
        b_zero[] <- 0
        if (distofempty == "zero" & any(rowSums(linknt.RewSha) == 0)) {
          b_os.raw <- b_zero
        }
        else {
          b_os.raw <- vegdist(linknt.RewSha, method = index, binary = binary)
        }
      }

      # calculate dissimilarity of species, (b_s), interactions (b_wn) and rewiring (b_os) (by betadiver)
      if (function.dist == "betadiver") {
        if (binary == FALSE) {
          warning("betadiver only uses binary data; for quantitative indices use vegdist")
        }
        else {
          b_s <- betadiver(specnt.all, method = index)
          b_wn <- betadiver(linknt, method = index)
          if (distofempty == "zero" & any(rowSums(linknt.RewSha) == 0)) {
            b_os.raw <- b_zero
          }
          else {
            b_os.raw <- betadiver(linknt.RewSha, method = index)
          }
        }
      }

      b_os <- b_os.raw
      # calculate species turnover component (b_st)
      b_st <- b_wn - b_os.raw
      return(c(WN = b_wn, S = b_s, OS = b_os, ST = b_st))
    }

    # execute "poisotMethod" for conventional approach
    pMethod.net <- poisotMethod(linknt.sharedsp.net, specnt.all.net)

    # if c1 and c2 exist, execute "poisotMethod for INIS approach
    if(!is.null(c1) && !is.null(c2)) {
      pMethod.com <- poisotMethod(linknt.sharedsp.com, specnt.all.com)
    }

    # Organize Poisot results
    output <- c("WN" = pMethod.net[[1]], "S.net" = pMethod.net[[2]], "OS.net" = pMethod.net[[3]], "ST.net" = pMethod.net[[4]])
    if(!is.null(c1) && !is.null(c2)) {
      output <- c(output, "S.com" = pMethod.com[[2]], "OS.com" = pMethod.com[[3]], "ST.com" = pMethod.com[[4]])
    }
  }

  # NOVOTNY METHOD (Novotny 2009, Frund 2021)
  if(partitioning == "commondenom") {

    commondenomMethod <- function(linknt.sharedsp, specnt.all) {
      # calculate a, b and c components for applying dissimilarity formula
      # A: number of shared interactions;
      # B.tot and C.tot: unique interactions to network 1 and 2 respectively;
      # B.rew and C.rew: unique interactions to network 1 and 2, respectively, involving shared species.
      # B.uni and C.uni: unique interactions to network 1 and 2, respectively, involving unique species.
      pmins <- pmin(linknt[1, ], linknt[2, ])
      A <- sum(pmins)
      B.tot <- sum(linknt[1, ] - pmins)
      B.rew <- sum(linknt.sharedsp[1, ] - pmin(linknt.sharedsp[1,], linknt.sharedsp[2, ]))
      B.uni <- B.tot - B.rew
      C.tot <- sum(linknt[2, ] - pmins)
      C.rew <- sum(linknt.sharedsp[2, ] - pmin(linknt.sharedsp[1,], linknt.sharedsp[2, ]))
      C.uni <- C.tot - C.rew

      # define the denominator of the formula according to the index used.
      if (index == "bray") {
        index <- "sorensen"
      }
      if (index == "sorensen") {
        denominator <- 2 * A + B.tot + C.tot
      }
      if (index == "jaccard") {
        denominator <- A + B.tot + C.tot
      }

      # calculate dissimilarity of species (b_s), interactions (b_wn) 
      # and rewiring (b_os) and species turnover (b_st)
      b_wn <- (B.tot + C.tot)/denominator
      b_os <- (B.rew + C.rew)/denominator
      b_st <- (B.uni + C.uni)/denominator
      b_s <- vegdist(specnt.all, method = switch(index, jaccard = "jaccard", sorensen = "bray"), binary = binary)

      # if partition.st = TRUE, estimate which part of b_st is due to changes in low trophic level (b_st.l), high trophic
      # level (b_st.h) or both (b_st.lh)
      if (partition.st) {

        nt1.sharedlow <- nt1
        nt2.sharedlow <- nt2

        for(row in 1:nrow(nt1.sharedlow)) {
          if((length(nt1[row,][nt1[row,] > 0]) > 0) != (length(nt2[row,][nt2[row,] > 0]) > 0)) {
            nt1.sharedlow[row,] <- rep(0, ncol(nt1.sharedlow))
            nt2.sharedlow[row,] <- rep(0, ncol(nt2.sharedlow))
          }
        }

        nt1.sharedhigh <- nt1
        nt2.sharedhigh <- nt2
        for(col in 1:ncol(nt1.sharedhigh)) {
          if((length(nt1[,col][nt1[,col] > 0]) > 0) != (length(nt2[,col][nt2[,col] > 0]) > 0)) {
            nt1.sharedhigh[,col] <- rep(0, nrow(nt1.sharedhigh))
            nt2.sharedhigh[,col] <- rep(0, nrow(nt2.sharedhigh))
          }
        }

        nt1.sharedhighORlow <- nt1
        nt2.sharedhighORlow <- nt2

        for(row in 1:nrow(nt1.sharedhighORlow)) {
          if((length(nt1[row,][nt1[row,] > 0]) > 0) != (length(nt2[row,][nt2[row,] > 0]) > 0)) {
            nt1.sharedhighORlow[row,] <- rep(0, ncol(nt1.sharedhighORlow))
            nt2.sharedhighORlow[row,] <- rep(0, ncol(nt2.sharedhighORlow))
          }
        }

        nt1.sharedhighORlow <- nt1.sharedlow
        for(row in 1:nrow(nt1.sharedhighORlow)) {
          for(col in 1:ncol(nt1.sharedhighORlow)) {
            if(nt1.sharedlow[row, col] == 0) {
              nt1.sharedhighORlow[row, col] <- nt1.sharedhigh[row, col]
            }
          }
        }

        nt2.sharedhighORlow <- nt2.sharedlow
        for(row in 1:nrow(nt2.sharedhighORlow)) {
          for(col in 1:ncol(nt2.sharedhighORlow)) {
            if(nt2.sharedlow[row, col] == 0) {
              nt2.sharedhighORlow[row, col] <- nt2.sharedhigh[row, col]
            }
          }
        }

        linknt.lh <- matrix2links(nt1 - nt1.sharedhighORlow, nt2 - nt2.sharedhighORlow)
        linknt.l <- matrix2links(nt1 - nt1.sharedlow, nt2 - nt2.sharedlow) - linknt.lh
        linknt.h <- matrix2links(nt1 - nt1.sharedhigh, nt2 - nt2.sharedhigh) - linknt.lh

        B.l <- sum(linknt.l[1, ])
        B.h <- sum(linknt.h[1, ])
        B.lh <- sum(linknt.lh[1, ])
        C.l <- sum(linknt.l[2, ])
        C.h <- sum(linknt.h[2, ])
        C.lh <- sum(linknt.lh[2, ])
        b_st.l <- (B.l + C.l)/denominator
        b_st.h <- (B.h + C.h)/denominator
        b_st.lh <- (B.lh + C.lh)/denominator
      }

      # if partition.rr = TRUE
      if (partition.rr) {
        if (proportions == TRUE) {
          warning("partitionining into replacement and richness (abundance) difference components may be meaningless with proportions")
        }
        b_wn.repl <- 2 * min(B.tot, C.tot)/denominator
        b_os.repl <- 2 * min(B.rew, C.rew)/denominator
        b_wn.rich <- abs(B.tot - C.tot)/denominator
        b_os.rich <- abs(B.rew - C.rew)/denominator
      }
      output <- c(WN = b_wn, S = b_s, OS = b_os,  ST = b_st)
      if (partition.st == TRUE) {
        output <- c(output, ST.l = b_st.l, ST.h = b_st.h,
                    ST.lh = b_st.lh)
      }
      if (partition.rr == TRUE) {
        output <- c(output, WN.repl = b_wn.repl, OS.repl = b_os.repl,
                    WN.rich = b_wn.rich, OS.rich = b_os.rich)
      }
      return(output)

    }

    # execute "commondenomMethod" for conventional approach
    cMethod.net <- commondenomMethod(linknt.sharedsp.net, specnt.all.net)

    # if c1 and c2 exist, execture "commondenomMethod" for INIS approach
    if(!is.null(c1) && !is.null(c2)) {
      cMethod.com <- commondenomMethod(linknt.sharedsp.com, specnt.all.com)
    }

    # Organize Novotny results
    output <- c("WN" = cMethod.net[[1]], "S.net" = cMethod.net[[2]], "OS.net" = cMethod.net[[3]], "ST.net" = cMethod.net[[4]])
    if(!is.null(c1) || !is.null(c2)) {
      if(is.null(c1) || is.null(c2)) {
        warning("Independent community data missing. Only returning network results")
      } else {
        output <- c(output, "S.com" = cMethod.com[[2]], "OS.com" = cMethod.com[[3]], "ST.com" = cMethod.com[[4]] )
      }
    }
    if(partition.rr) {
      output <- c(output, "WN.repl" = cMethod.net[[8]], "WN.rich" = cMethod.net[[10]], "OS.net.repl" = cMethod.net[[9]], "OS.net.rich" = cMethod.net[[11]])
      if(!is.null(c1) && !is.null(c2)) {
        output <- c(output, "OS.com.repl" = cMethod.com[[9]], "OS.com.rich" = cMethod.com[[11]])
      }
    }
    if(partition.st) {
      output <- c(output, "ST.net.l" = cMethod.net[[5]], "ST.net.h" = cMethod.net[[6]], "ST.net.lh" = cMethod.net[[7]])
      if(!is.null(c1) && !is.null(c2)) {
        output <- c(output, "ST.com.l" = cMethod.com[[5]], "ST.com.h" = cMethod.com[[6]], "ST.com.lh" = cMethod.com[[7]])
      }
    }
  }
  
  return(output)
}

# END OF SCRIPT
