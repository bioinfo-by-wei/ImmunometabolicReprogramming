library(dplyr)
library(stringi)
library(stringr)
library(Matrix)
library(osqp)

gene_num <- METAFlux:::gene_num
Hgem <- METAFlux:::Hgem
iso <- METAFlux:::iso
multi_comp <- METAFlux:::multi_comp
simple_comp <- METAFlux:::simple_comp

parallel_calculate_avg_exp <- function(myseurat,myident,n_bootstrap,seed) {
  set.seed(seed)
  samples=METAFlux:::generate_boots(myseurat@meta.data[,myident],n_bootstrap)
  exp <- mclapply(1:n_bootstrap,METAFlux:::get_ave_exp,myseurat,samples,myident)
  exp <- do.call(cbind, exp)
  return(exp)
}

# parallel_calculate_avg_exp <- function(myseurat,myident,n_bootstrap,seed) {
#   set.seed(seed)
#   samples=METAFlux:::generate_boots(myseurat@meta.data[,myident],n_bootstrap)
#   cl <- makeCluster(24, type = "PSOCK")
#   exp <- parLapply(cl,1:n_bootstrap,METAFlux:::get_ave_exp,myseurat,samples,myident)
#   stopCluster(cl)
#   exp <- do.call(cbind, exp)
#   return(exp)
# }

stdize = function(x, ...) {(x ) / (max(x, ...))}

parallel_calculate_reaction_score <- function (data){
    if (sum(data < 0) > 0)
        stop("Expression data needs to be all positive")
    features <- rownames(data)
    if (sum(features %in% rownames(gene_num)) == 0)
        stop("Requested gene names cannot be found. Rownames of input data should be human gene names.Please check the rownames of input data.")
    message(paste0(round(sum(features %in% rownames(gene_num))/3625 *
        100, 3), "% metabolic related genes were found......"))
    message("Computing metabolic reaction activity scores......")
    core <- do.call(rbind, mclapply(1:length(iso), METAFlux:::calculate_iso_score,
        data = data, list = iso, gene_num = gene_num))
    core2 <- do.call(rbind, mclapply(1:length(simple_comp), METAFlux:::calculate_simple_comeplex_score,
        data, list = simple_comp, gene_num = gene_num))
    core3 <- do.call(rbind, mclapply(1:length(multi_comp), METAFlux:::calculate_multi_comp,
        data = data, gene_num = gene_num))
    message("Preparing for score matrix......")
    rownames(core) <- names(iso)
    rownames(core2) <- names(simple_comp)
    rownames(core3) <- names(multi_comp)
    big_score_matrix <- rbind(core, core2, core3)
    big_score_matrix <- apply(big_score_matrix, 2, stdize, na.rm = T)
    big_score_matrix[is.na(big_score_matrix)] <- 0
    empty_helper <- as.data.frame(Hgem$Reaction)
    colnames(empty_helper) <- "reaction"
    Final_df <- merge(empty_helper, big_score_matrix, all.x = T,
        by.x = 1, by.y = 0)
    Final_df[is.na(Final_df)] <- 1
    rownames(Final_df) <- Final_df$reaction
    Final_df$reaction <- NULL
    Final_df <- Final_df[Hgem$Reaction, , drop = F]
    if (all.equal(rownames(Final_df), Hgem$Reaction)) {
        message("Metabolic reaction activity scores successfully calculated \n")
    }
    else {
        message("Calculation not reliable Check input data format \n")
    }
    Final_df[which(Hgem$LB == 0 & Hgem$UB == 0), ] <- 0
    return(Final_df)
}

fast_compute_sc_flux <- function (num_cell, fraction, fluxscore, medium){
    if (!all.equal(sum(fraction),1))
        stop("Sum of fractions must be eqaul to 1")
    if (length(fraction) != num_cell)
        stop("Number of cell clusters does not match with the length of fraction")
    Hgem <- METAFlux:::Hgem
    mat <- Hgem$S
    reaction_name <- Hgem$Reaction
    names(reaction_name) <- NULL
    A_combined <- METAFlux:::A_combined
    D <- Matrix(0,8378,13082,sparse=T)
    candi_list <- lapply(rapply(list(mat, lapply(1:2, function(x) D)),
        enquote, how = "unlist"), eval)
    matrix_construct <- diag(1, ncol = num_cell, nrow = num_cell)
    matrix_construct[matrix_construct == 0] <- 2
    celltype_matrix <- NULL
    A_matrix <- NULL
    for (i in 1:num_cell) {
        A_matrix[[i]] <- A_combined
        celltype_matrix[[i]] <- do.call(cbind, candi_list[matrix_construct[,
            i]])
    }
    message("Preparing for TME S matrix.....")
    final_s <- rbind(do.call(cbind, A_matrix), do.call(rbind,
        celltype_matrix))
    whole3 <- rbind(as.sparse(diag(-1, 1648, 1648)), as.sparse(matrix(0,
        num_cell * nrow(mat), 1648)))
    final_s <- cbind(final_s, whole3)
    external <- paste("external_medium", reaction_name[which(Hgem$pathway ==
        "Exchange/demand reactions")])
    reaction_name[which(Hgem$pathway == "Exchange/demand reactions")] <- paste0("internal_medium ",
        reaction_name[which(Hgem$pathway == "Exchange/demand reactions")])
    cell_reaction <- NULL
    for (i in 1:num_cell) {
        cell_reaction[[i]] <- paste(paste("celltype", i), reaction_name)
    }
    construct_reaction_names <- c(unlist(cell_reaction), external)
    P1 <- as.sparse(Diagonal(ncol(final_s),1))
    message("S matrix completed......")
    flux_vector <- list()
    message("Compute metabolic flux......")
    Seq <- seq(1, ncol(fluxscore), by = num_cell)
    pb <- txtProgressBar(0, length(Seq), style = 3)
    for (t in Seq) {
        setTxtProgressBar(pb, match(t, Seq))
        message("  ",t,' in ',max(Seq))
        score <- fluxscore[, c(t:(t + num_cell - 1))]
        #P <- as.sparse(diag(c(unlist(mcMap(rep, fraction, 13082)),rep(1, 1648)), ncol(final_s), ncol(final_s)))
        P <- as.sparse(Diagonal(ncol(final_s),c(unlist(Map(rep, fraction, 13082)),rep(1, 1648))))
        fraction_finals <- final_s %*% P
        q <- rep(0, ncol(final_s))
        q[seq(13015, ncol(final_s), by = 13082)] <- -10000 *
            fraction
        A <- as.sparse(rbind(fraction_finals, P1))
        ras <- c(as.vector(unlist(score[, 1:num_cell])), rep(1,
            1648))
        origlb <- c(rep(Hgem$LB, num_cell), rep(-1, 1648))
        origlb[rep(Hgem$rev == 1, num_cell)] <- (-ras[rep(Hgem$rev ==
            1, num_cell)])
        origlb[rep(Hgem$rev == 0, num_cell)] <- 0
        origub <- ras
        origlb[tail(1:ncol(final_s), 1648)] <- 0
        matches <- sapply(medium$reaction_name, function(x) intersect(which(stri_detect_fixed(construct_reaction_names,
            x)), tail(1:ncol(final_s), 1648)))
        origlb[matches] <- -1
        l <- c(rep(0, nrow(final_s)), origlb)
        u <- c(rep(0, nrow(final_s)), origub)
        settings <- osqpSettings(max_iter = 1000000L, eps_abs = 1e-04,linsys_solver = 1,
            eps_rel = 1e-04, verbose = FALSE, adaptive_rho_interval = 50)
        model <- osqp(P, q, A, l, u, settings)
        res <- model$Solve()
        flux_vector <- append(flux_vector, list(res$x))
    }
    close(pb)
    flux_vector <- as.data.frame(do.call(cbind, flux_vector))
    rownames(flux_vector) <- construct_reaction_names
    return(flux_vector)
}

