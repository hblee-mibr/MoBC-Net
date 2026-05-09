# #' module distance function
# #' 
# #' Calculate the distance and z_score of the communities
# #' use.f can be choosed in the DistFunction.r
# #' ex) get.shortest.dist, get.kernel.dist, get.centre.dist, get.separation.dist, get.closest.dist
# #' As a default, get.closest.dist function is used to measure the distance between communities
# #' User also can make the dist function and use it for calculating module distance
# #' This function is made to know the z-score of a measured distance from distances of degree-preserved random networks
# #' 


# # randomization 1
# simple_sampling <- function(deg, membership, ncv = 1000) {

#     tnodes <- 1:length(membership)
#     m = sum(membership==1)
#     n = sum(membership==2)
#     sind <- sample(tnodes, m + n, replace = FALSE)

#     Urand <- sind[1:m]
#     Wrand <- sind[(m + 1):(m + n)]
# 	return(list(Urand, Wrand))
# }





# #----- randomization 2

# abinning_estimate_deg_bag_prob <- function(deg, membership, ncv = 1000, kappa = 0.07) {
#     stopifnot(length(deg) == length(membership))
#     # m: 모듈 개수(최대 인덱스), 질병 모듈: membership > 0
#     m <- max(membership, na.rm = TRUE)
#     disease_module_idx <- which(membership > 0)

#     # 고유 degree(오름차순)별로 초기 bin 구성
#     xtick <- sort(unique(deg))
#     L <- length(xtick)

#     xbag <- vector("list", L)        # 각 bin이 포함하는 degree 값
#     node_bag <- vector("list", L)    # 각 bin이 포함하는 노드 인덱스
#     kk <- matrix(0L, nrow = m, ncol = L)  # 모듈별 카운트(질병 모듈만 집계)
#     nn <- rep(Inf, L)                # 각 bin의 max(min_k)*ncv  (min_k = tkk/tnn)

#     # tabulate helper: membership 값 1..m 에 대해 카운트
#     count_by_module <- function(mod_vec, m) {
#     # mod_vec은 1..m 범위여야 함 (membership>0만 전달)
#     if (length(mod_vec) == 0) return(integer(m))
#     tab <- tabulate(mod_vec, nbins = m)
#     as.integer(tab)
#     }

#     # 초기 bin별 통계 계산
#     for (i in seq_along(xtick)) {
#     xbag[[i]] <- xtick[i]
#     node_bag[[i]] <- which(deg == xtick[i])

#     # 질병 모듈에 속한 노드만 선택
#     tsel_node_bag <- intersect(node_bag[[i]], disease_module_idx)

#     # membership 1..m 에 대한 히스토그램
#     tkk <- count_by_module(membership[tsel_node_bag], m)
#     kk[, i] <- tkk

#     mind <- which(tkk > 0L)
#     if (length(mind) > 0L) {
#         tnn <- length(node_bag[[i]])           # 전체 노드 수(비질병 포함)
#         min_k <- tkk[mind] / tnn               # 각 모듈 비율
#         nn[i] <- max(min_k) * ncv              # 해당 bin의 "최대 모듈 비율 × ncv"
#     }
#     }

#     threshold <- ncv * kappa

#     # 병합 루프: 모든 nn < threshold가 될 때까지 오른쪽에서부터 위반 bin 병합
#     while (!all(nn < threshold)) {
#     viol <- which(nn >= threshold)
#     pt <- max(viol)   # 가장 오른쪽 위반 bin

#     if (pt != 1L) {
#         # 왼쪽 bin(pt-1)과 병합
#         xbag[[pt - 1L]] <- c(xbag[[pt - 1L]], xbag[[pt]])
#         node_bag[[pt - 1L]] <- c(node_bag[[pt - 1L]], node_bag[[pt]])
#         kk[, pt - 1L] <- kk[, pt - 1L] + kk[, pt]

#         # nn 재계산 (질병 모듈이 하나라도 있으면 갱신)
#         tkk <- kk[, pt - 1L]
#         mind <- which(tkk > 0L)
#         if (length(mind) > 0L) {
#         tnn <- length(node_bag[[pt - 1L]])
#         min_k <- tkk[mind] / tnn
#         nn[pt - 1L] <- max(min_k) * ncv
#         }
#     } else {
#         # 첫 번째 bin이면 오른쪽(pt+1)과 병합
#         xbag[[pt + 1L]] <- c(xbag[[pt]], xbag[[pt + 1L]])
#         node_bag[[pt + 1L]] <- c(node_bag[[pt]], node_bag[[pt + 1L]])
#         kk[, pt + 1L] <- kk[, pt] + kk[, pt + 1L]

#         tkk <- kk[, pt + 1L]
#         mind <- which(tkk > 0L)
#         if (length(mind) > 0L) {
#         tnn <- length(node_bag[[pt + 1L]])
#         min_k <- tkk[mind] / tnn
#         nn[pt + 1L] <- max(min_k) * ncv
#         }
#     }

#     # 병합된 pt 제거 (리스트/행렬/벡터 동기화)
#     keep <- setdiff(seq_along(xbag), pt)
#     xbag <- xbag[keep]
#     node_bag <- node_bag[keep]
#     kk <- kk[, keep, drop = FALSE]
#     nn <- nn[keep]
#     }

#     n <- vapply(node_bag, length, integer(1))
#     k <- kk
#     #   list(xbag = xbag, node_bag = node_bag, n = n, k = k)
#     list(deg_bag = xbag, node_bag = node_bag, n = n, k = k)
# }



# # estimate_deg_bag <- function(deg, membership, ratiov, ncv = 1000) {
# #   # ncv보다 작게 들어오면 기본값으로 설정
# #   # (MATLAB: if nargin < 3, ncv = 1000;)
# #   # R에서는 기본 인자로 처리하므로 별도 처리 생략
  
# #   # threshold 설정
# #   threshold <- ncv * ratiov
  
# #   # 고유한 deg 값 (MATLAB의 unique(deg))
# #   xtick <- unique(deg) %>% sort # unique degree
  
# #   # 리스트(셀 배열) 초기화
# #   deg_bag <- vector("list", length(xtick))
# #   node_bag <- vector("list", length(xtick))
# #   sel_node_bag <- vector("list", length(xtick))
# #   ncases <- numeric(length(xtick))
  
# #   # deg_bag, node_bag, sel_node_bag 구성
# #   for (i in seq_along(xtick)) { # index
# #     deg_bag[[i]] <- xtick[i] #묶여진 degree들
# #     node_bag[[i]] <- which(deg == xtick[i]) # degree가 같은 node index 담아
# #     sel_node_bag[[i]] <- which(deg == xtick[i] & membership > 0)  # degree가 같으면서 u,w cluster에 속한 node index 담아
    
# #     Ni <- length(node_bag[[i]])
# #     ni <- length(sel_node_bag[[i]])
    
# #     # ncases(i) = nchoosek(Ni, ni)
# #     # R에서는 choose 함수를 사용
# #     if (Ni >= ni) {
# #       ncases[i] <- choose(Ni, ni)
# #     } else {
# #       # 예외적으로 Ni < ni이면 0이나 NA 처리 등 필요
# #       ncases[i] <- 0
# #     }
# #   }
  
# #   # index 설정
# #   # - ind_noselect: sel_node_bag[[i]]가 비어있는(선택할 노드가 없는) 위치
# #   # - ind_select: sel_node_bag[[i]]가 비어있지 않은(선택해야 할 노드가 있는) 위치
# #   ind_noselect <- which(sapply(sel_node_bag, length) == 0)
# #   ind_select <- which(sapply(sel_node_bag, length) > 0)
  
# # # (ind_select의 마지막 + 1)부터 끝까지는 noselect 대상에서 제외
# #   if (length(ind_select) > 0) {
# #     last_sel <- tail(ind_select, 1)
# #     # MATLAB의 setdiff(ind_noselect, (ind_select(end)+1):length(sel_node_bag))
# #     remove_range <- seq(last_sel + 1, length(sel_node_bag))
# #     ind_noselect <- setdiff(ind_noselect, remove_range)
# #   }
  
# #   # noselect인 곳은 ncases를 ncv로 설정
# #   ncases[ind_noselect] <- ncv
  
# #   # 첫 번째 while 루프
# #   # MATLAB: while prod(ncases > threshold) ~= 1
# #   # R: while (!all(ncases > threshold))
# #   while (!all(ncases >= threshold)) {
# #     # pt = max(find(ncases < threshold)) : 마지막으로 조건 미충족인 인덱스
# #     pt <- max(which(ncases < threshold))
# #     # cat(pt,'\n')
# #     # pt가 1이 아닌 경우 이전 bag과 병합, 1인 경우 다음 bag과 병합
# #     if (pt != 1) {
# #       deg_bag[[pt - 1]] <- c(deg_bag[[pt - 1]], deg_bag[[pt]])
# #       node_bag[[pt - 1]] <- c(node_bag[[pt - 1]], node_bag[[pt]])
# #       sel_node_bag[[pt - 1]] <- c(sel_node_bag[[pt - 1]], sel_node_bag[[pt]])
# #     } else {
# #       # pt가 1인 경우
# #       if (pt + 1 <= length(deg_bag)) {
# #         deg_bag[[pt + 1]] <- c(deg_bag[[pt]], deg_bag[[pt + 1]])
# #         node_bag[[pt + 1]] <- c(node_bag[[pt]], node_bag[[pt + 1]])
# #         sel_node_bag[[pt + 1]] <- c(sel_node_bag[[pt]], sel_node_bag[[pt + 1]])
# #       }
# #     }
    
# #     # 병합한 pt번째 bag 제거
# #     keep_idx <- setdiff(seq_along(deg_bag), pt)
# #     deg_bag <- deg_bag[keep_idx]
# #     node_bag <- node_bag[keep_idx]
# #     sel_node_bag <- sel_node_bag[keep_idx]
    
# #     # ncases 재계산
# #     ncases <- mapply(function(nb, sb) {
# #       Ni <- length(nb)
# #       ni <- length(sb)
# #       if (Ni >= ni) choose(Ni, ni) else 0
# #     }, node_bag, sel_node_bag)

# # # noselect / select 재설정
# #     ind_noselect <- which(sapply(sel_node_bag, length) == 0)
# #     ind_select <- which(sapply(sel_node_bag, length) > 0)
# #     if (length(ind_select) > 0) {
# #       last_sel <- tail(ind_select, 1)
# #       remove_range <- seq(last_sel + 1, length(sel_node_bag))
# #       ind_noselect <- setdiff(ind_noselect, remove_range)
# #     }
# #     ncases[ind_noselect] <- ncv
    
# #     # MATLAB: display(num2str(length(node_bag)))
# #     # cat("현재 node_bag의 길이:", length(node_bag), "\n")
# #   }
  
# #   # 두 번째 while 루프
# #   # node_bag 내 원소들의 길이가 단조감소(diff_nn <= 0)여야 함
# #   nn <- sapply(node_bag, length)
# #   diff_nn <- diff(nn)
  
# #   # MATLAB: while prod(diff_nn <= 0) ~= 1
# #   while (!all(diff_nn <= 0)) {
# #     # pt = max(find(diff_nn>0)) : 단조감소가 깨지는 마지막 인덱스
# #     pt <- max(which(diff_nn > 0))
    
# #     # pt가 1이 아닌 경우 이전 bag과 병합
# #     if (pt != 1) {
# #       deg_bag[[pt - 1]] <- c(deg_bag[[pt - 1]], deg_bag[[pt]])
# #       node_bag[[pt - 1]] <- c(node_bag[[pt - 1]], node_bag[[pt]])
# #       sel_node_bag[[pt - 1]] <- c(sel_node_bag[[pt - 1]], sel_node_bag[[pt]])
# #     } else {
# #       # pt가 1인 경우 다음 bag과 병합
# #       if (pt + 1 <= length(deg_bag)) {
# #         deg_bag[[pt + 1]] <- c(deg_bag[[pt]], deg_bag[[pt + 1]])
# #         node_bag[[pt + 1]] <- c(node_bag[[pt]], node_bag[[pt + 1]])
# #         sel_node_bag[[pt + 1]] <- c(sel_node_bag[[pt]], sel_node_bag[[pt + 1]])
# #       }
# #     }
    
# #     # pt번째 bag 제거
# #     keep_idx <- setdiff(seq_along(deg_bag), pt)
# #     deg_bag <- deg_bag[keep_idx]
# #     node_bag <- node_bag[keep_idx]
# #     sel_node_bag <- sel_node_bag[keep_idx]
    
# #     # diff 재계산
# #     nn <- sapply(node_bag, length)
# #     diff_nn <- diff(nn)
    
# #     # cat("현재 deg_bag의 길이:", length(deg_bag), "\n")
# #   }
# #     cat("final deg_bag의 길이:", length(deg_bag), "\n")
  
# #   # 결과 반환
# #   return(list(deg_bag = deg_bag, node_bag = node_bag))
# # }




# # randomization 3
# clustermn <- function(S, m, n) {
#   #-----------------------------------------------------
#   # S: p x p 거리 행렬 (대각 0, 대칭)
#   # m: 첫 번째 클러스터의 노드 개수
#   # n: 두 번째 클러스터의 노드 개수 (m + n = p)
#   #-----------------------------------------------------
  
# 	p <- nrow(S)

# 	if (m > p || m < 1) {
# 	stop("m은 1과 p 사이여야 합니다.")
# 	}
# 	if (m + n != p) {
# 	stop("m + n이 행렬 S의 크기 p와 같아야 합니다.")
# 	}

# 	# 1. 거리 행렬 S를 가우시안 커널로 변환 -> 유사도 행렬 S2
# 	# sigma = 거리 행렬 S의 전체 평균
# 	sigma <- mean(S)
# 	# exp(-S^2 / sigma^2) 연산: 행렬 원소별로 계산
# 	# R에서 ^ 연산은 원소별이므로, S^2는 원소별 제곱
# 	S2 <- exp(-(S^2) / (sigma^2))

# 	# 2. 그래프 라플라시안 L = D - S2
# 	# D는 각 행 합(=노드 차수)을 대각에 배치한 행렬
# 	deg <- rowSums(S2)
# 	D <- diag(deg)
# 	L <- D - S2

# 	# 3. 라플라시안의 고유값 분해 (Fiedler vector: 두 번째로 작은 고유값에 대응하는 고유벡터)
# 	# R의 eigen()은 (대체로) 고유값을 내림차순으로 정렬하여 반환하므로,
# 	# 오름차순 정렬 인덱스를 사용해서 두 번째 작은 고유값에 해당하는 고유벡터를 찾는다.
# 	ev <- eigen(L, symmetric = TRUE)
# 	idx_asc <- order(ev$values)              # 고유값 오름차순 정렬 인덱스
# 	fiedler <- ev$vectors[, idx_asc[2]]      # 두 번째로 작은 고유값의 고유벡터

# 	# 4. Fiedler vector 값 기준 오름차순 정렬
# 	sortedIdx <- order(fiedler)

# 	# 5. 앞에서부터 m개 노드를 cluster 1, 나머지를 cluster 2로 할당
# 	clusterAssignment <- integer(p)  # p개 길이의 정수 벡터(초기 0)
# 	clusterAssignment[sortedIdx[1:m]] <- 1
# 	clusterAssignment[sortedIdx[(m + 1):p]] <- 2

# 	# 필요한 경우, 결과 확인 용 출력 가능:
# 	# print(clusterAssignment)

# 	return(clusterAssignment)
# }






# # Randomization 4. 

# modularity_sampling<-function(re, deg, membership, S){
# 	cnode_bag = re$node_bag
# 	names(cnode_bag) = 1:length(cnode_bag)
# 	cdeg_bag = re$deg_bag


# 	### (a) U 모듈 무작위 생성
# 	Udeg <- deg[membership==1]               # U 내 노드들의 차수
# 	Udeg <- sort(Udeg, decreasing = TRUE)
		
# 	Urand <- c()                    # 최종 뽑힌 U 모듈 노드
		
# 	for (i in seq_along(Udeg)) {
# 		# cat(i,'\n')
# 		# xbag 중 현재 차수를 포함하는 인덱스를 찾는다.
# 		ibag <- which(sapply(cdeg_bag, function(x) any(x %in% Udeg[i])))

# 		if (length(Urand) == 0) {
# 			# 첫 노드 샘플링
# 			tUr <- sample(cnode_bag[[ibag]], 1, replace = FALSE)
# 			# 뽑힌 노드는 제거
# 			cnode_bag[[ibag]] <- setdiff(cnode_bag[[ibag]], tUr)
# 			# Ur 갱신
# 			Urand <- c(Urand, tUr)
# 		} else {
# 			# Urand 안에 이미 노드가 있다면, 그 노드들과의 거리에 따라 가중치 설정
# 			# S[Urand, cnode_bag[[ibag]]] 부분행렬의 열방향 평균(colMeans)이 각 후보 노드까지의 평균거리
# 			tdist <- colMeans(S[Urand, cnode_bag[[ibag]], drop = FALSE])
			
# 			# 가우시안 형태의 가중치 (tsigma=1 고정)
# 			tsigma <- 1
# 			weights <- exp(-(tdist^2) / (2 * (tsigma^2)))
# 			weights <- weights / sum(weights)
			
# 			# 가중치를 이용해 노드 하나 샘플링
#             if(length(weights)==1){
#         		tUr <- cnode_bag[[ibag]] #sample(cnode_bag[[ibag]], 1, replace = FALSE, prob = weights)
#             } else{
#     			tUr <- sample(cnode_bag[[ibag]], 1, replace = FALSE, prob = weights)
#             }
# 			# tUr <- sample(cnode_bag[[ibag]], 1, replace = FALSE, prob = weights)
# 			# 뽑힌 노드는 제거
# 			cnode_bag[[ibag]] <- setdiff(cnode_bag[[ibag]], tUr)
# 			# Ur 갱신
# 			Urand <- c(Urand, tUr)
# 		}
# 		}


#     ### (b) W 모듈 무작위 생성
#     Wdeg <- deg[membership==2]               # U 내 노드들의 차수
#     Wdeg <- sort(Wdeg, decreasing = TRUE)
    
#     Wrand <- c()
    
#     for (i in seq_along(Wdeg)) {
# 		ibag <- which(sapply(cdeg_bag, function(x) any(x %in% Wdeg[i])))

# 		if (length(Wrand) == 0) {
# 			tWr <- sample(cnode_bag[[ibag]], 1, replace = FALSE)
# 			cnode_bag[[ibag]] <- setdiff(cnode_bag[[ibag]], tWr)
# 			Wrand <- c(Wrand, tWr)
# 		} else {
# 			# 이미 Wr가 있을 경우엔 가우시안 가중치로 샘플링
# 			tdist <- colMeans(S[Wrand, cnode_bag[[ibag]], drop = FALSE])
			
# 			tsigma <- 1
# 			weights <- exp(-(tdist^2) / (2 * (tsigma^2)))
# 			weights <- weights / sum(weights)

#             if(length(weights)==1){
#         		tWr <- cnode_bag[[ibag]] #sample(cnode_bag[[ibag]], 1, replace = FALSE, prob = weights)        
#             } else{
#     			tWr <- sample(cnode_bag[[ibag]], 1, replace = FALSE, prob = weights)
#             }
# 			cnode_bag[[ibag]] <- setdiff(cnode_bag[[ibag]], tWr)
# 			Wrand <- c(Wrand, tWr)
# 		}
#     }
# 	return(list(Urand, Wrand))
# }


# modularity_sampling_multi <- function(re, deg, membership, S) {
#     cnode_bag = re$node_bag
#     names(cnode_bag) = 1:length(cnode_bag)
#     cdeg_bag = re$deg_bag
    
#     # 모든 모듈 ID 가져오기
#     module_ids <- unique(membership[membership > 0]) %>% sort
#     result_list <- list()
    
#     for (mod_id in module_ids) {
#         # 현재 모듈의 차수들
#         mod_deg <- deg[membership == mod_id]
#         # mod_deg <- sort(mod_deg, decreasing = TRUE)
#         # cat(mod_id,'\n')
#         # cat('hh')
#         mod_deg <- sample(mod_deg, length(mod_deg), replace=FALSE)
#         # cat('kk\n')

#         mod_rand <- c()
#         for (i in seq_along(mod_deg)) {
#             ibag <- which(sapply(cdeg_bag, function(x) any(x %in% mod_deg[i])))
            
#             if (length(mod_rand) == 0) {
#                 # 첫 노드 샘플링
#                 selected_node <- sample(cnode_bag[[ibag]], 1, replace = FALSE)
#                 cnode_bag[[ibag]] <- setdiff(cnode_bag[[ibag]], selected_node)
#                 mod_rand <- c(mod_rand, selected_node)
#             } else {
#                 # 가중치 기반 샘플링
#                 tdist <- colMeans(S[mod_rand, cnode_bag[[ibag]], drop = FALSE])
#                 weights <- exp(-(tdist^2) / (2 * 1^2))
#                 weights <- weights / sum(weights)
                
#                 if(length(weights)==1){
#                     selected_node = cnode_bag[[ibag]]
#                 } else{
#                     selected_node <- sample(cnode_bag[[ibag]], 1, replace = FALSE, prob = weights)
#                 }
#                 cnode_bag[[ibag]] <- setdiff(cnode_bag[[ibag]], selected_node)
#                 mod_rand <- c(mod_rand, selected_node)
#             }
#         }
        
#         result_list[[mod_id]] <- mod_rand
#     }
    
#     return(result_list)
# }


