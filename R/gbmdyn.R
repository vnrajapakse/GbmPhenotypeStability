map_sample_exp_to_pathway_states <- function(exp_mat, pathways) {
  set.seed(1)
  sample_pwstates <- GSVA::gsva(expr = exp_mat, gset.idx.list = pathways, 
                                method = 'ssgsea', verbose=FALSE) %>% t() %>% as_tibble()
  return(sample_pwstates)
}

match_samples_to_attrs <- function(sample_pwstates, attr_pwstates, 
                                   sample_ids=NULL, pwattr_ids = NULL) {
  #set.seed(1)
  stopifnot(identical(colnames(sample_pwstates), colnames(attr_pwstates)))
  if (!is.null(sample_ids)) {
    stopifnot(length(sample_ids) == nrow(sample_pwstates))
  }
  if (!is.null(pwattr_ids)) {
    stopifnot(length(pwattr_ids) == nrow(attr_pwstates))
  }
  
  sample_attr_tbl <- tibble(
    sample_uid = seq_len(nrow(sample_pwstates)),
    pwattr_uid = NA,
    t_on_vs_off = NA,
    pval = NA
  )
  
  tbl_row_to_vec <- function(k, tb) {
    set_names(as.numeric(tb[k, ]), colnames(tb))
  }
  
  for (i in seq_len(nrow(sample_pwstates))) {
    sample_i <- tbl_row_to_vec(i, sample_pwstates)
    attr_match_pvals  <- rep(NA, nrow(attr_pwstates))
    attr_match_tstats <- rep(NA, nrow(attr_pwstates))
    
    for (j in seq_len(nrow(attr_pwstates))) {
      attr_j <- tbl_row_to_vec(j, attr_pwstates)
      
      tt_result <- tryCatch(
        t.test(sample_i[attr_j == 1], sample_i[attr_j == 0], alternative = 'greater'), 
        error = function(e) { return(NA) }
      )
      
      if (is.list(tt_result)) {
        attr_match_pvals[j]  <- tt_result$p.value
        attr_match_tstats[j] <- tt_result$statistic 
      }
    }
    
    j_match <- which.min(attr_match_pvals)
    if (!is.na(j_match)) {
      sample_attr_tbl[[i, 'pwattr_uid']] <- j_match
      sample_attr_tbl[[i, 't_on_vs_off']] <- attr_match_tstats[[j_match]]
      sample_attr_tbl[[i, 'pval']] <- attr_match_pvals[[j_match]]
    }
  }
  
  if (!is.null(sample_ids)) {
    sample_attr_tbl[['sample_uid']] <- sample_ids[sample_attr_tbl[['sample_uid']]]
  }
  if (!is.null(pwattr_ids)) {
    sample_attr_tbl[['pwattr_uid']] <- pwattr_ids[sample_attr_tbl[['pwattr_uid']]]
  }
  
  return(sample_attr_tbl)
}


vnorm <- function(x) {
  sqrt(sum(x^2))
}