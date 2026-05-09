


#--------------------------- new


get.freq <-function(g, snode, enode){
	edges = igraph::all_shortest_paths(g, snode, enode)
	# edges = edges$res %>% lapply(function(xx) setdiff(names(xx), names(xx)[c(1,length(xx))]))
    # edges = edges$res %>% lapply(function(xx) names(xx[2:(length(xx)-1)]))
    # edges = edges$res %>% lapply(function(xx) names(xx[2:(length(xx)-1)]))
    edges = edges$res %>% lapply(function(xx) if (length(xx) > 2) names(xx[2:(length(xx)-1)]) else character(0))
	# etab = edges %>% unlist
	return(edges) #!!
}

get.freq.v2 <-function(g, snode, enode){
	edges = igraph::all_shortest_paths(g, snode, enode)
	# edges = edges$res %>% lapply(function(xx)  setdiff(names(xx), names(xx)[1]))
	# edges = edges$res %>% lapply(function(xx)  names(xx[2:length(xx)]))
	edges = edges$res %>% lapply(function(xx)  if (length(xx) > 2) names(xx[2:(length(xx)-1)]) else character(0))
	# etab = edges %>% unlist
	return(edges) #!!
}


#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title Get.Centrality
#' @param network results from CommDistFunction function
#' @param module1 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @param module2 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @returns data.frame
#' @export
#' @examples
#' Get.Centrality(MoBC.result, 'module_1','module_2')

MoBC.genes <- function(network,
                             module.gene.list,
                             module1, module2,
                            #  randomMethod=c('None','RandC','RandCD','RandCM','RandCDM'),
                            randomMethod=c('None','RandSD','RandSDM'),
							 random = 1000,
                             nCore=1,
                             ratio = 0.1) {
    overlap_filtering=TRUE
    # cat(method,'\n')
    if (is.character(randomMethod)){
        randomMethod <- match.arg(randomMethod)
        # cat(dist.function,'\n')
        randomMethod <- switch(randomMethod,
            None = 'None',
            # RandC = 'random1',
            RandSD = 'randSD',
            # RandCM = 'random3',
            RandSDM = 'RandSDM')
    } else {
        stop('Method function is wrong. Check the method function', call.=FALSE)
    }
	if (!any(names(module.gene.list) %in% module1) | !any(names(module.gene.list) %in% module2)) {
		stop('Module name is not included in module list. Please assign right module names.', call.=FALSE)
	}
	if (any(lengths(module.gene.list)==0)) {
		stop('Please assign right node names', call.=FALSE)
	}
	if (!is.character(module1) | !is.character(module2)) {
		stop('Please assign character name as module name', call.=FALSE)
	}

	g.res  <- preprocessedNetwork(network)
    comm.genelist <- CommunityGenelist(module.gene.list, g.res, overlap_filtering = overlap_filtering)


	x=cal.MoBCgenes(g.res, comm.genelist,
					community1n=module1, 
					community2n=module2,
                    random=random, ratio=ratio,randomMethod=randomMethod, nCore=nCore)
    # x = subset(x, score>0)
	return(x)
	}




#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title cal.MoBCgenes
#' @param g graph
#' @param comm.genelist list of community genes
#' @param community1n The name of the community for which centrality is being calculated. This should be one of the communities provided as input
#' @param community2n The name of the community for which centrality is being calculated. This should be one of the communities provided as input
#' @returns vector
#' @export
#' @examples
#' cal.MoBCgenes(graph, 'community1','community2',random,ratio,randomMethod, nCore)



cal.MoBCgenes <- function(g, comm.genelist, community1n, community2n,random,ratio,randomMethod, nCore){
# >
    # cat('cal.MoBCgenes')

    community1 = comm.genelist[[community1n]]
    community2 = comm.genelist[[community2n]]

    allg = igraph::V(g)$name %>% as.character()
    scorev = cal.MoBCgenes.values(g, community1, community2, allg)

	score.df = data.frame(gene=allg,score=scorev)
	score.df$node_type = 'link'
	score.df$node_type[score.df$gene %in% c(community1, community2)] = 'community genes'
	score.df = score.df %>% dplyr::arrange(-score)
    score.df = subset(score.df, score>0)

    colix = c('gene','score','node_type')
    # cat('test.\n')
    if(randomMethod!='None'){
        cat('Normalized MoBC score will be provided.\n')
        colix = c('gene','normalized_score','node_type')
        random.mat = cal.MoBC.random(g, comm.genelist, community1n, community2n,random,ratio,randomMethod,show.binning=FALSE, nCore=nCore)
        pval = sapply(score.df$gene, function(gn){
            xval = score.df[match(gn, score.df$gene),'score']
            pval = sum(random.mat[gn,]>=xval)/(random+1)
        })
        nscore = sapply(score.df$gene, function(gn){
            xval = score.df[match(gn, score.df$gene),'score']
            zvalv = c(random.mat[gn,],xval)
            (xval-mean(zvalv))/sd(zvalv)
        })
        score.df$normalized_score = nscore
        score.df$pval = pval
        # pv = sapply(score.df$pval, function(vv) ifelse(vv>0.5, 1-vv,vv))
        score.df$p.adj = p.adjust(pval,'BH')
        colix = c('gene','score','normalized_score','node_type','pval','p.adj')

    }
	return(score.df[,colix])
}



#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title cal.MoBCgenes.values
#' @param g graph
#' @param community1 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @param community2 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @returns vector
#' @export
#' @examples
#' cal.MoBCgenes.values(graph, 'module_1','module_2', allg)


cal.MoBCgenes.values <- function(g, community1, community2, allg){
    
	scorevec = rep(0, length(igraph::V(g))) %>% 'names<-'(igraph::V(g)$name)
	shortestm = igraph::distances(g, community1, community2)
	rmin  = apply(shortestm,1,function(xx) colnames(shortestm)[which(xx %in% min(xx))])
	cmin  = apply(shortestm,2,function(xx) rownames(shortestm)[which(xx %in% min(xx))])


    r.endNode = lengths(rmin) %>% sum
    c.endNode = lengths(cmin) %>% sum

    # comm1
    r.sp.genel = sapply(names(rmin), function(start.node){

		end.node = rmin[[start.node]]
        edges = igraph::all_shortest_paths(g, start.node, end.node)
        end.ix = edges$res %>% sapply(function(xx) names(xx)[length(xx)])
        end.ix = split(1:length(end.ix), end.ix)

        vv = sapply(end.ix, function(ixix){
            # intg = edges$res[ixix] %>% lapply(function(xx) setdiff(names(xx), names(xx)[c(1,length(xx))]))
            etab = edges$res[ixix] %>% lapply(function(xx) if (length(xx) > 2) names(xx[2:(length(xx)-1)]) else character(0))
            etab1 = etab %>% unlist %>% table
            val = etab1 / length(etab)
            if (length(etab1) == 0) {
                val = setNames(rep(0, length(allg)), allg)
                # val[names(etab1)] = as.numeric(etab1) / length(etab)
            }
    		return(val[allg])
        })%>% 'rownames<-'(allg)
        vv[is.na(vv)] = 0

        # vv = sapply(end.node, function(endn){
    	# 	etab = get.freq(g, start.node, endn)
        #     etab1 = etab %>% unlist %>% table
        #     val = etab1/length(etab)
    	# 	return(val[allg])
        # }) %>% 'rownames<-'(allg)
        # vv[is.na(vv)] = 0
        return(vv)
	})
    rval = do.call(cbind, r.sp.genel)
    if(ncol(rval)!=r.endNode) stop('Not matched column number (r)', call. = FALSE)
    r.val = apply(rval,1,sum)

    # comm2
    c.sp.genel = sapply(names(cmin), function(start.node){

		end.node = cmin[[start.node]]
        edges = igraph::all_shortest_paths(g, start.node, end.node)
        end.ix = edges$res %>% sapply(function(xx) names(xx)[length(xx)])
        end.ix = split(1:length(end.ix), end.ix)

        vv = sapply(end.ix, function(ixix){
            # intg = edges$res[ixix] %>% lapply(function(xx) setdiff(names(xx), names(xx)[c(1,length(xx))]))
            etab = edges$res[ixix] %>% lapply(function(xx) if (length(xx) > 2) names(xx[2:(length(xx)-1)]) else character(0))
            etab1 = etab %>% unlist %>% table
            val = etab1 / length(etab)
            if (length(etab1) == 0) {
                val = setNames(rep(0, length(allg)), allg)
                # val[names(etab1)] = as.numeric(etab1) / length(etab)
            }
    		return(val[allg])
        })%>% 'rownames<-'(allg)
        vv[is.na(vv)] = 0

        # vv = sapply(end.node, function(endn){
    	# 	etab = get.freq(g, start.node, endn)
        #     etab1 = etab %>% unlist %>% table
        #     val = etab1/length(etab)
    	# 	return(val[allg])
        # }) %>% 'rownames<-'(allg)
        # vv[is.na(vv)] = 0
        return(vv)
	})
    cval = do.call(cbind, c.sp.genel)
    if(ncol(cval)!=c.endNode) stop('Not matched column number (r)', call. = FALSE)
    c.val = apply(cval,1,sum)

    scorev = (r.val+c.val)/sum(r.endNode+c.endNode) #!!

	# score.df = data.frame(gene=allg,community1.score =as.numeric(r.val), community2.score=as.numeric(c.val), score=all.score)
    # scorev = score.df$community1.score + score.df$community2.score
    # names(scorev) = allg

	return(scorev)
}

# #---->> old,,, so takes too long
# cal.MoBCgenes.values <- function(g, community1, community2, allg){
    
# 	scorevec = rep(0, length(igraph::V(g))) %>% 'names<-'(igraph::V(g)$name)
# 	shortestm = igraph::distances(g, community1, community2)
# 	rmin  = apply(shortestm,1,function(xx) colnames(shortestm)[which(xx %in% min(xx))])
# 	cmin  = apply(shortestm,2,function(xx) rownames(shortestm)[which(xx %in% min(xx))])


#     r.endNode = lengths(rmin) %>% sum
#     c.endNode = lengths(cmin) %>% sum

#     # comm1
#     r.sp.genel = sapply(names(rmin), function(start.node){
# 		end.node = rmin[[start.node]]
#         vv = sapply(end.node, function(endn){
#     		etab = get.freq(g, start.node, endn)
#             etab1 = etab %>% unlist %>% table
#             val = etab1/length(etab)
#     		return(val[allg])
#         }) %>% 'rownames<-'(allg)
#         vv[is.na(vv)] = 0
#         return(vv)
# 	})
#     rval = do.call(cbind, r.sp.genel)
#     if(ncol(rval)!=r.endNode) stop('Not matched column number (r)', call. = FALSE)
#     r.val = apply(rval,1,sum)

#     # comm2
#     c.sp.genel = sapply(names(cmin), function(start.node){
# 		end.node = cmin[[start.node]]
#         vv = sapply(end.node, function(endn){
#     		etab = get.freq(g, start.node, endn)
#             etab1 = etab %>% unlist %>% table
#             val = etab1/length(etab)
#     		return(val[allg])
#         }) %>% 'rownames<-'(allg)
#         vv[is.na(vv)] = 0
#         return(vv)
# 	})
#     cval = do.call(cbind, c.sp.genel)
#     if(ncol(cval)!=c.endNode) stop('Not matched column number (r)', call. = FALSE)
#     c.val = apply(cval,1,sum)



#     scorev = (r.val+c.val)/sum(r.endNode+c.endNode) #!!

# 	# score.df = data.frame(gene=allg,community1.score =as.numeric(r.val), community2.score=as.numeric(c.val), score=all.score)
#     # scorev = score.df$community1.score + score.df$community2.score
#     # names(scorev) = allg

# 	return(scorev)
# }




# g = res2@graph
# community1 = res2@filtered.communities[[2]]
# community2 = res2@filtered.communities[[3]]
# random = 1000
# ratio = 0.1

# cal.MoBC.random(g, comm.genelist, community1n, community2n,random,ratio,cal.p,show.binning=FALSE, nCore=nCore)
# community1 --> genes in community1
# community2 --> genes in community2
# comm.genelist


# random part
cal.MoBC.random <- function(g, comm.genelist, community1n, community2n,random,ratio,randomMethod,show.binning, nCore){
    # cat("cal.MoBC.radnom - T_T")
    allg = igraph::V(g)$name %>% as.character()

    cl1g = comm.genelist[[community1n]]
    cl2g = comm.genelist[[community2n]]


    fn = make_cache_key(g, modules=comm.genelist, params=paste0(c(random,ratio),collapse='_'))
    dirn = paste0('./MoBCtmp/',fn,'/',randomMethod)

    # ㅡ make membership
    deg = igraph::degree(g)
    membership <- rep(0, length(deg))
    for(ii in 1:length(comm.genelist)){
        membership[names(deg) %in% comm.genelist[[ii]]] <- ii #membership --> numeric
    }

    # - make file name
    module_names <- names(comm.genelist)
    files_valid <- check_module_files(dirn, module_names, random)
    cat(files_valid,'\n')
    cat(randomMethod,'\n')


    # make directory
    if(!all(files_valid)){
        cat(paste0("You don't have optimal tmp files for random sampling - ",randomMethod,". Generating random samples will take time...  \n"))

        if(file.exists(dirn)) unlink(dirn, recursive=TRUE)
        dir.create(dirn, recursive=TRUE,showWarnings = FALSE)
        hist.bin0 = abinning_estimate_deg_bag_prob(deg, membership, kappa=ratio, ncv = random)

        fn1 = paste0(dirn,'/hist_bin.RDS')
        saveRDS(hist.bin0, file=fn1)
    }


    # if not random --> make random file first
    if(!all(files_valid) & randomMethod=='randSD'){
        hist.bin = hist.bin0$node_bag
        names(hist.bin) = 1:length(hist.bin)
        make_random = parallel::mclapply(1:random, mc.cores=nCore, function(j){

            sample.save = list()
            for(ii in 1:length(comm.genelist)){
                # cat('Random samples for ',names(comm.genelist)[ii],' module\n')
                samplingN = sapply(hist.bin, function(xx) sum(names(xx) %in% comm.genelist[[ii]])) %>% 'names<-'(names(hist.bin))
                clg.random = lapply(names(samplingN), function(xn){
                    use.bg = setdiff(hist.bin[[xn]], unlist(sample.save)) 
                    sample(use.bg, samplingN[[xn]], replace=FALSE)
                }) %>% unlist %>% unique

                sample.save[[ii]] = clg.random
                dir.create(paste0(dirn,'/',names(comm.genelist)[ii]), recursive=TRUE,showWarnings = FALSE)
                write.csv(clg.random,paste0(dirn,'/',names(comm.genelist)[ii],'/rand',j,'.csv'), row.names=F)
            }
            return(NULL)
        })

    } else if(!all(files_valid) & randomMethod=='RandSDM'){
        S <- igraph::distances(g, algorithm = "unweighted")

        comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
            rsamplel = modularity_sampling_multi(hist.bin0, deg, membership, S)     
            for(ii in 1:length(rsamplel)){
                # cat('Random samples for ',names(comm.genelist)[ii],' module\n')

                clg.random = rsamplel[[ii]]
                dir.create(paste0(dirn,'/',names(comm.genelist)[ii]), recursive=TRUE,showWarnings = FALSE)
                write.csv(clg.random,paste0(dirn,'/',names(comm.genelist)[ii],'/rand',j,'.csv'), row.names=F)
            }
            return(NULL)
        })

    } else if(randomMethod!='RandSD' & randomMethod!='RandSDM'){
        stop('Please enter the right method for randomization.', call.=FALSE)
    }


    # check again
    files_valid1 <- check_module_files(dirn, module_names, random)

    if(all(files_valid1)){### edit needed

        cat(paste0('You have tmp files for random sampling - ',randomMethod,". We will use these files.\n"))

        comm.distance.list = parallel::mclapply(1:random,mc.cores=nCore, function(j){
        # comm.distance.list = lapply(1:random,function(j){
            cat('random ',j,'\n')
            fnn=paste0(dirn,'/',community1n,'_',community2n,'/MoBC_rand',j,'_ix.RDS')
            if(!file.exists(paste0(dirn,'/',community1n,'_',community2n))) dir.create(paste0(dirn,'/',community1n,'_',community2n))
            if(file.exists(fnn)){
                comm.distance = readRDS(file=fnn)
                return(comm.distance)
            }
            rs1 = read.csv(paste0(dirn,'/',community1n,'/rand',j,'.csv'))[,1]
            rs2 = read.csv(paste0(dirn,'/',community2n,'/rand',j,'.csv'))[,1]
            
            comm.distance = cal.MoBCgenes.values(g, rs1,rs2, allg) 
            saveRDS(comm.distance,file=fnn)
            return(comm.distance)
        })
        comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)
        return(comm.distance.list)
    }

    cat('fin comm.dist\n')
    # start.time-end.time
    return(comm.distance.list)
}




#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title FCS.genes
#' @param network results from CommDistFunction function
#' @param module1 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @param module2 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @returns data.frame
#' @export
#' @examples
#' FCS.genes(MoBC.result, 'module_1','module_2')

FCS.genes <- function(network,
                             module.gene.list,
                             module1, module2,
                            #  randomMethod=c('None','RandC','RandCD','RandCM','RandCDM'),
                            randomMethod=c('None','RandSD','RandSDM'),
							 random = 1000,
                             nCore=1,
                             ratio = 0.1) {
    overlap_filtering=TRUE
    # cat(method,'\n')
    if (is.character(randomMethod)){
        randomMethod <- match.arg(randomMethod)
        # cat(dist.function,'\n')
        randomMethod <- switch(randomMethod,
            None = 'None',
            # RandC = 'random1',
            RandSD = 'randSD',
            # RandCM = 'random3',
            RandSDM = 'RandSDM')
    } else {
        stop('Method function is wrong. Check the method function', call.=FALSE)
    }
	if (!any(names(module.gene.list) %in% module1) | !any(names(module.gene.list) %in% module2)) {
		stop('Module name is not included in module list. Please assign right module names.', call.=FALSE)
	}
	if (any(lengths(module.gene.list)==0)) {
		stop('Please assign right node names', call.=FALSE)
	}
	if (!is.character(module1) | !is.character(module2)) {
		stop('Please assign character name as module name', call.=FALSE)
	}

	g.res  <- preprocessedNetwork(network)
    comm.genelist <- CommunityGenelist(module.gene.list, g.res, overlap_filtering = overlap_filtering)

    cat("We are running FCS\n")
	x=cal.FCSgenes(g.res, comm.genelist,
					community1n=module1, 
					community2n=module2,
                    random=random, ratio=ratio,randomMethod=randomMethod, nCore=nCore)
    # x = subset(x, score>0)
	return(x)
	}





# random part
cal.random <- function(g, comm.genelist, community1n, community2n,random,ratio,randomMethod,show.binning, nCore, method='FCS'){
    # cat("cal.FCS.radnom - T_T")
    allg = igraph::V(g)$name %>% as.character()

    cl1g = comm.genelist[[community1n]]
    cl2g = comm.genelist[[community2n]]


    fn = make_cache_key(g, modules=comm.genelist, params=paste0(c(random,ratio),collapse='_'))
    dirn = paste0('./MoBCtmp/',fn,'/',randomMethod)

    # ㅡ make membership
    deg = igraph::degree(g)
    membership <- rep(0, length(deg))
    for(ii in 1:length(comm.genelist)){
        membership[names(deg) %in% comm.genelist[[ii]]] <- ii #membership --> numeric
    }

    # - make file name
    module_names <- names(comm.genelist)
    files_valid <- check_module_files(dirn, module_names, random)
    cat(files_valid,'\n')
    cat(randomMethod,'\n')


    # make directory
    if(!all(files_valid)){
        cat(paste0("You don't have optimal tmp files for random sampling - ",randomMethod,". Generating random samples will take time...  \n"))

        if(file.exists(dirn)) unlink(dirn, recursive=TRUE)
        dir.create(dirn, recursive=TRUE,showWarnings = FALSE)
        hist.bin0 = abinning_estimate_deg_bag_prob(deg, membership, kappa=ratio, ncv = random)

        fn1 = paste0(dirn,'/hist_bin.RDS')
        saveRDS(hist.bin0, file=fn1)
    }


    # if not random --> make random file first
    if(!all(files_valid) & randomMethod=='randSD'){
        hist.bin = hist.bin0$node_bag
        names(hist.bin) = 1:length(hist.bin)
        make_random = parallel::mclapply(1:random, mc.cores=nCore, function(j){

            sample.save = list()
            for(ii in 1:length(comm.genelist)){
                # cat('Random samples for ',names(comm.genelist)[ii],' module\n')
                samplingN = sapply(hist.bin, function(xx) sum(names(xx) %in% comm.genelist[[ii]])) %>% 'names<-'(names(hist.bin))
                clg.random = lapply(names(samplingN), function(xn){
                    use.bg = setdiff(hist.bin[[xn]], unlist(sample.save)) 
                    sample(use.bg, samplingN[[xn]], replace=FALSE)
                }) %>% unlist %>% unique

                sample.save[[ii]] = clg.random
                dir.create(paste0(dirn,'/',names(comm.genelist)[ii]), recursive=TRUE,showWarnings = FALSE)
                write.csv(clg.random,paste0(dirn,'/',names(comm.genelist)[ii],'/rand',j,'.csv'), row.names=F)
            }
            return(NULL)
        })

    } else if(!all(files_valid) & randomMethod=='RandSDM'){
        S <- igraph::distances(g, algorithm = "unweighted")

        comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
            rsamplel = modularity_sampling_multi(hist.bin0, deg, membership, S)     
            for(ii in 1:length(rsamplel)){
                # cat('Random samples for ',names(comm.genelist)[ii],' module\n')

                clg.random = rsamplel[[ii]]
                dir.create(paste0(dirn,'/',names(comm.genelist)[ii]), recursive=TRUE,showWarnings = FALSE)
                write.csv(clg.random,paste0(dirn,'/',names(comm.genelist)[ii],'/rand',j,'.csv'), row.names=F)
            }
            return(NULL)
        })

    } else if(randomMethod!='RandSD' & randomMethod!='RandSDM'){
        stop('Please enter the right method for randomization.', call.=FALSE)
    }


    # check again
    files_valid1 <- check_module_files(dirn, module_names, random)

    if(all(files_valid1)){### edit needed

        cat(paste0('You have tmp files for random sampling - ',randomMethod,". We will use these files.\n"))

        if(method=='FCS'){
            cat('You choosed FCS method.\n')
            comm.distance.list = parallel::mclapply(1:random,mc.cores=nCore, function(j){
                cat(j,'-th random sampling is done.\n')
                fnn=paste0(dirn,'/',community1n,'_',community2n,'/FCS_rand',j,'_ix.RDS')
                if(file.exists(fnn)){
                    comm.distance = readRDS(file=fnn)
                    return(comm.distance)
                }
                # comm.distance.list = lapply(1:random,function(j){
                rs1 = read.csv(paste0(dirn,'/',community1n,'/rand',j,'.csv'))[,1]
                rs2 = read.csv(paste0(dirn,'/',community2n,'/rand',j,'.csv'))[,1]
                
                comm.distance = cal.FCSgenes.values(g, rs1,rs2, allg) 
                saveRDS(comm.distance,file=fnn)
                return(comm.distance)
            })
            comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)
            return(comm.distance.list)

        } else if(method=='S2B'){

            cat('You choosed S2B method.\n')
            S <- igraph::distances(g, algorithm = "unweighted")
            avgd <- igraph::mean_distance(g)
            Smt = S<avgd
            if(!file.exists(paste0(dirn,'/',community1n,'_',community2n))) dir.create(paste0(dirn,'/',community1n,'_',community2n))
            for(j in 1:random){
                fnn=paste0(dirn,'/',community1n,'_',community2n,'/S2B_rand',j,'_ix.RDS')
                if(file.exists(fnn)) next
                cat(j,'-th random path is saved....\n')
                rs1 = read.csv(paste0(dirn,'/',community1n,'/rand',j,'.csv'))[,1] #%>% as.character
                rs2 = read.csv(paste0(dirn,'/',community2n,'/rand',j,'.csv'))[,1] #%>% as.character
                use_ids = apply(Smt[rs1,],1,function(xx) intersect(which(xx), rs2))
                saveRDS(use_ids,file=fnn)
            }
            rm(S)
            rm(Smt)

            comm.distance.list = parallel::mclapply(1:random,mc.cores=nCore, function(j){
            # comm.distance.list = lapply(1:random,function(j){
                cat(j,'-th random sampling is done.\n')
                use_ids = readRDS(file=paste0(dirn,'/',community1n,'_',community2n,'/S2B_rand',j,'_ix.RDS'))
                rs1 = read.csv(paste0(dirn,'/',community1n,'/rand',j,'.csv'))[,1] #%>% as.character
                rs2 = read.csv(paste0(dirn,'/',community2n,'/rand',j,'.csv'))[,1] #%>% as.character
                comm.distance = cal.S2Bgenes.values(g, rs1,rs2, allg,use_ids=use_ids)  #!!!!!!!!!!!!!
                return(comm.distance)
            })
            comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)
            return(comm.distance.list)

        }

    }

    cat('fin comm.dist\n')
    # start.time-end.time
    return(comm.distance.list)
}





#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title cal.FCSgenes.values
#' @param g graph
#' @param community1 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @param community2 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @returns data.frame
#' @export
#' @examples
#' cal.FCgene(graph, 'module_1','module_2')


cal.FCSgenes.values <- function(g, community1, community2, allg){

    gene.ix = igraph::V(g)$name

    # re = igraph::all_shortest_paths(g, community1[1:3],community2[1:4])
    # lapply(re$res, function(xx) names(xx)[1] ) %>% unlist %>% unique
    # lapply(re$res, function(xx) names(xx)[length(xx)] ) %>% unlist %>% unique
    
    re = sapply(community1, function(g1){
            
        edges = igraph::all_shortest_paths(g, g1, community2)
        end.ix = edges$res %>% sapply(function(xx) names(xx)[length(xx)])
        end.ix = split(1:length(end.ix), end.ix)
        resl = lapply(end.ix, function(ixix){
            # intg = edges$res[ixix] %>% lapply(function(xx) setdiff(names(xx), names(xx)[c(1,length(xx))]))
            intg = edges$res[ixix] %>% lapply(function(xx) if (length(xx) > 2) names(xx[2:(length(xx)-1)]) else character(0))
            intg.tab = unlist(intg) %>% table
            intg.tab = intg.tab/length(intg)
            intg.tab[gene.ix] %>% 'names<-'(gene.ix)
        })
        res = do.call(rbind, resl) %>% as.matrix
        res[is.na(res)] = 0
        resv = apply(res,2,sum)
    })
    re1 = apply(re,1,sum)
    re1 = re1/length(community1)/length(community2)

	return(re1[allg])
}






#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title cal.FCSgenes
#' @param g graph
#' @param community1n The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @param community2n The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @returns data.frame
#' @export
#' @examples
#' cal.FCSgenes(g, comm.genelist, community1n, community2n,random,ratio,randomMethod, nCore)


cal.FCSgenes <- function(g, comm.genelist, community1n, community2n,random,ratio,randomMethod, nCore){

    # cat('cal.MoBCgenes')

    community1 = comm.genelist[[community1n]]
    community2 = comm.genelist[[community2n]]

    allg = igraph::V(g)$name %>% as.character()
    scorev = cal.FCSgenes.values(g, community1, community2, allg)

	score.df = data.frame(gene=allg,score=scorev)
	score.df$node_type = 'link'
	score.df$node_type[score.df$gene %in% c(community1, community2)] = 'community genes'
	score.df = score.df %>% dplyr::arrange(-score)
    score.df = subset(score.df, score>0)

    colix = c('gene','score','node_type')
    # cat('test.\n')
    if(randomMethod!='None'){
        cat('Normalized FCS score will be provided.\n')
        colix = c('gene','normalized_score','node_type')
        random.mat = cal.random(g, comm.genelist, community1n, community2n,random,ratio,randomMethod,show.binning=FALSE, nCore=nCore,method='FCS')
        pval = sapply(score.df$gene, function(gn){
            xval = score.df[match(gn, score.df$gene),'score']
            pval = sum(random.mat[gn,]>=xval)/(random+1)
        })
        nscore = sapply(score.df$gene, function(gn){
            xval = score.df[match(gn, score.df$gene),'score']
            zvalv = c(random.mat[gn,],xval)
            (xval-mean(zvalv))/sd(zvalv)
        })
        score.df$normalized_score = nscore
        score.df$pval = pval
        # pv = sapply(score.df$pval, function(vv) ifelse(vv>0.5, 1-vv,vv))
        score.df$p.adj = p.adjust(pval,'BH')
        colix = c('gene','score','normalized_score','node_type','pval','p.adj')

    }
	return(score.df[,colix])
}





#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title S2B.genes
#' @param network results from CommDistFunction function
#' @param module1 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @param module2 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @returns data.frame
#' @export
#' @examples
#' S2B.genes(MoBC.result, 'module_1','module_2')

S2B.genes <- function(network,
                             module.gene.list,
                             module1, module2,
                            #  randomMethod=c('None','RandC','RandCD','RandCM','RandCDM'),
                            randomMethod=c('None','RandSD','RandSDM'),
							 random = 1000,
                             nCore=1,
                             ratio = 0.1) {
    overlap_filtering=TRUE
    # cat(method,'\n')
    if (is.character(randomMethod)){
        randomMethod <- match.arg(randomMethod)
        # cat(dist.function,'\n')
        randomMethod <- switch(randomMethod,
            None = 'None',
            # RandC = 'random1',
            RandSD = 'randSD',
            # RandCM = 'random3',
            RandSDM = 'RandSDM')
    } else {
        stop('Method function is wrong. Check the method function', call.=FALSE)
    }
	if (!any(names(module.gene.list) %in% module1) | !any(names(module.gene.list) %in% module2)) {
		stop('Module name is not included in module list. Please assign right module names.', call.=FALSE)
	}
	if (any(lengths(module.gene.list)==0)) {
		stop('Please assign right node names', call.=FALSE)
	}
	if (!is.character(module1) | !is.character(module2)) {
		stop('Please assign character name as module name', call.=FALSE)
	}

	g.res  <- preprocessedNetwork(network)
    comm.genelist <- CommunityGenelist(module.gene.list, g.res, overlap_filtering = overlap_filtering)

    cat("We are running S2B\n")
	x=cal.S2Bgenes(g.res, comm.genelist,
					community1n=module1, 
					community2n=module2,
                    random=random, ratio=ratio,randomMethod=randomMethod, nCore=nCore)
    # x = subset(x, score>0)
	return(x)
	}



#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title cal.S2Bgenes.values
#' @param g graph
#' @param community1 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @param community2 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @returns data.frame
#' @export
#' @examples
#' cal.FCgene(graph, 'module_1','module_2')


cal.S2Bgenes.values <- function(g, community1, community2, allg, use_ids=NULL){
    
    gene.ix = igraph::V(g)$name
    if(is.null(use_ids)){
        S = igraph::distances(g, v = community1, to = community2, algorithm = "unweighted")
        avgd <- igraph::mean_distance(g)
        Smt = S<avgd
    
        re = sapply(community1, function(g1){
            # cat(g1,'\n')
            use.comm2 = community2[Smt[g1,community2]]
            if(length(use.comm2)==0) return(rep(0, length(allg)))
            edges = igraph::all_shortest_paths(g, g1, use.comm2)
            end.ix = edges$res %>% sapply(function(xx) names(xx)[length(xx)])
            end.ix = split(1:length(end.ix), end.ix)
            resl = lapply(end.ix, function(ixix){
                intg = edges$res[ixix] %>% lapply(function(xx) if (length(xx) > 2) names(xx[2:(length(xx)-1)]) else character(0))
                intg.tab = unlist(intg) %>% unique
                valval = rep(1, length(intg.tab)) %>% 'names<-'(intg.tab)
                valval[gene.ix] %>% 'names<-'(gene.ix)
            })
            res = do.call(rbind, resl) %>% as.matrix
            res[is.na(res)] = 0
            # six = Smt[g1,community2]
            resv = apply(as.data.frame(res),2,sum)
        })
        pool = sum(Smt)
    } else{
        re = sapply(1:length(community1), function(g1.ix){
            g1 = community1[[g1.ix]]
            # cat(g1,'\n')
            use.comm2 = use_ids[[g1.ix]]
            if(length(use.comm2)==0) return(rep(0, length(allg)))
            edges = igraph::all_shortest_paths(g, g1, use.comm2)
            end.ix = edges$res %>% sapply(function(xx) names(xx)[length(xx)])
            end.ix = split(1:length(end.ix), end.ix)
            resl = lapply(end.ix, function(ixix){
                intg = edges$res[ixix] %>% lapply(function(xx) if (length(xx) > 2) names(xx[2:(length(xx)-1)]) else character(0))
                intg.tab = unlist(intg) %>% unique
                valval = rep(1, length(intg.tab)) %>% 'names<-'(intg.tab)
                valval[gene.ix] %>% 'names<-'(gene.ix)
            })
            res = do.call(rbind, resl) %>% as.matrix
            res[is.na(res)] = 0
            # six = Smt[g1,community2]
            resv = apply(as.data.frame(res),2,sum)
        })
        pool = use_ids %>% unlist %>% length
    }
    re1 = apply(re,1,sum)

    if(max(re1)>pool) cat('ERRORO!!!!!!')
	return(re1/pool)
}




#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title cal.S2Bgenes
#' @param g graph
#' @param community1n The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @param community2n The name of the module for which centrality is being calculated. This should be one of the communities provided as input
#' @returns data.frame
#' @export
#' @examples
#' cal.S2Bgenes(g, comm.genelist, community1n, community2n,random,ratio,randomMethod, nCore)


cal.S2Bgenes <- function(g, comm.genelist, community1n, community2n,random,ratio,randomMethod, nCore){

    # cat('cal.MoBCgenes')

    community1 = comm.genelist[[community1n]]
    community2 = comm.genelist[[community2n]]

    allg = igraph::V(g)$name %>% as.character()
    scorev = cal.S2Bgenes.values(g, community1, community2, allg)

	score.df = data.frame(gene=allg,score=scorev)
	score.df$node_type = 'link'
	score.df$node_type[score.df$gene %in% c(community1, community2)] = 'community genes'
	score.df = score.df %>% dplyr::arrange(-score)
    score.df = subset(score.df, score>0)

    colix = c('gene','score','node_type')
    # cat('test.\n')
    if(randomMethod!='None'){
        cat('Normalized S2B score will be provided.\n')
        colix = c('gene','normalized_score','node_type')
        random.mat = cal.random(g, comm.genelist, community1n, community2n,random,ratio,randomMethod,show.binning=FALSE, nCore=nCore,method='S2B')
        pval = sapply(score.df$gene, function(gn){
            xval = score.df[match(gn, score.df$gene),'score']
            pval = sum(random.mat[gn,]>=xval)/(random+1)
        })
        nscore = sapply(score.df$gene, function(gn){
            xval = score.df[match(gn, score.df$gene),'score']
            zvalv = c(random.mat[gn,],xval)
            (xval-mean(zvalv))/sd(zvalv)
        })
        score.df$normalized_score = nscore
        score.df$pval = pval
        # pv = sapply(score.df$pval, function(vv) ifelse(vv>0.5, 1-vv,vv))
        score.df$p.adj = p.adjust(pval,'BH')
        colix = c('gene','score','normalized_score','node_type','pval','p.adj')

    }
	return(score.df[,colix])
}



#' Calculate centrality between two modules from MoBC result 
#' 
#' 
#' @title plotDist
#' @param MoBC.result results from CommDistFunction function
#' @param pval cut-off for filtering edges between communities
#' @returns plot
#' @export
#' @examples
#' plotDist(MoBC.result, pval=0.05)




plotDist <- function(MoBC.result, pval=0.05){
	if(!is(MoBC.result, 'MoBCresult')){
		stop("input should be MoBC class", call. = FALSE)
	}

	distm = MoBC.result@MoBCresults
	sig.dist = subset(distm, pvalue < pval)[,1:3]
	sig.dist$weight = -sig.dist$z_score
	ntkg = igraph::graph_from_data_frame(sig.dist[,c('Module1','Module2','weight')], directed=FALSE)
	ntkg = igraph::simplify(ntkg, remove.multiple = TRUE, remove.loops = TRUE)

    maxn = max(lengths(MoBC.result@filtered.modules))
    comm.col = colorspace::sequential_hcl(length(MoBC.result@filtered.modules), "Terrain") %>% 'names<-'(names(MoBC.result@filtered.modules))
    
    commn = lengths(MoBC.result@filtered.modules)[igraph::V(ntkg)$name]
    sizev = (commn-min(commn))/(max(commn)-min(commn))
    sizev = sizev*20+20

	layout <- igraph::layout_with_fr(ntkg)

	plre = plot(ntkg, 
		layout = layout, 
		# mark.groups = split(V(g)$name,clv),
		# vertex.label = fgid1[match(V(g)$name, fgid1$EntrezID),'gene_name'],
		# vertex.label = '', #vns
		vertex.color=comm.col[igraph::V(ntkg)$name],
		vertex.frame.width=0.3,
		vertex.frame.color='white',
		edge.color ="grey",#adjustcolor('black', alpha=0.6),
		# vertex.size= (cln[V(cl.ntkg)$name]^0.5)*4,
		vertex.size=sizev,
		# vertex.label.dist=1,
		# vertex.frame.color = 'grey90',
		vertex.label.color='black',
		# vertex.label.font=ifelse(V(g)$name %in% np.gl[[pn]], 2,1),
		vertex.label.size = 0.1,
		edge.width=(rank(igraph::E(ntkg)$weight))
	)

    legend("bottomright", col=comm.col, pch=19, legend=names(comm.col), title='Module')
}



