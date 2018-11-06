find_min_projection <- function(x, number) {
	return(min(unlist(lapply(x@curves, function(y) {
		return(y$dist[number]) 
	} ))))
}

res <- list()
middle_element <- list()
sampling <- FALSE
#Get original clusters number without sampling
source("clusterCurves.r")
clusters_orig <- ipmcMD@ident
sampling <- TRUE
names_t <- colnames(ipmcMD@data)
names(clusters_orig) <- gsub("X", "",names(clusters_orig))
#Here we start sampling
for (i in seq(1, 100)) {
tryCatch ({
	#Get new clusters number
	source("clusterCurves.r")
	clusters <- ipmcMD@ident
	new_cluster_names<-names(clusters)
	#Collect minimal distances between curves and cells
	for (i in seq_along(new_cluster_names)) {
		if (is.null(middle_element[new_cluster_names[i]])) {
			middle_element[[new_cluster_names[i]]] <- list()
		}
		middle_element[[new_cluster_names[i]]] <- c(middle_element[[new_cluster_names[i]]], find_min_projection(slingObjMD, i))
	}
	#Collect all clusters number in this round
	full_t <- lapply(seq_along(names_t), function(j) {
		x <- names_t[j]
		if (x %in% names(clusters)) {
			print(as.character(clusters[[x]]))
			return(as.character(clusters[[x]]))
		}
		else {
			return("0")
		}
	})
	names(full_t) <- names_t#gsub(" ", ".", names_t)
	sub_table <- list()
	for (i in 1:(max(as.numeric(full_t))+1)) {
		a<-names(which(full_t == toString(i)))
		max_intersect <- 0
		max_b <- 0
		for (j in 1:(max(as.numeric(clusters_orig))+1)) {
			b<-names(clusters_orig)[which(clusters_orig == toString(j))]
			if (max_intersect < length(intersect(a,b))) {
				max_intersect <- length(intersect(a,b))
				max_b <- j
			}
		}
		sub_table <- c(sub_table, max_b)
	}
	#If we don't have cluster in this round, we mark cluster number as zero
	new_full_t <- lapply(full_t, function(x) {
		if (as.numeric(x) == 0) {
			return("0")
		}
		else {
			return(sub_table[[as.numeric(x)]])
	}})
	#Write to file
	fileConn<-file("table_permutation.txt", 'a')
	sub_line <- ""
	for (i in 3:length(sub_table)-1) {
		if (i == 2) {
			sub_line<-paste(sub_line,sub_table[[i]], sep="")
		}
		else {
			sub_line<-paste(sub_line,sub_table[[i]], sep=" ")
		}
	}
	writeLines(sub_line, fileConn)
	close(fileConn)
	names(new_full_t) <- names(full_t)
	res <- c(res,  new_full_t)
}, error = function(e) {
	print(e);
	print("It's sad, but we will continue")})
}
names_uniq <- sort(unique(names(res)))
names_res <- unique(names(res))
res_final <- list()
for (i in seq(1, length(names_res))) {
	print(names_res[i])
	res_final[[names_res[i]]] <- as.numeric(unlist(res[which(names(res) == names_res[i])]))
}
sink("output.txt")
for (i in seq(1, length(names(res_final)))) {
	cat(c(names(res_final)[[i]],res_final[[i]]))
	cat('\n')
}
sink()
#Find mean and sd for distances between curve and cell
sd_m<-unlist(lapply(middle_element, function(x) {return(sd(x))}))
mean_m<-unlist(lapply(middle_element, function(x) {return(mean(x))}))
res <- data.frame(MEAN=mean_m, SD=sd_m)
q()
