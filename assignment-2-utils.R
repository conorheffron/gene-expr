import_data <- function(filename, pattern_regex, skip_n) {
  data_files = list.files(paste(filename, "/", sep = ""), pattern = pattern_regex)
  data_sets <- list()
  for (file in data_files) {
    if (endsWith(file, '.txt')) {
      print(paste(file, "- importing data"))
      data_sets[[file_path_sans_ext(file)]] <-
        import(paste(filename, file, sep = "/"),
               skip = skip_n,
               format = "txt")
    } else {
      print(paste(file, "is not needed for import..."))
    }
  }
  return(data_sets)
}

count_agg <-
  function(df,
           grp_col,
           n_results = 20,
           digits = 0,
           remove_ind = T) {
    if (is.null(df) | is.null(grp_col)) {
      return(NULL)
    }
    col_name <- ensym(grp_col)
    grp_df <- df |> group_by(!!col_name) |> count()
    grp_df_prop <-
      grp_df |> mutate(Freq = round(100 * n / sum(grp_df$n), digits))
    grp_o <- order(grp_df_prop[["n"]], decreasing = TRUE)
    grp <- data.frame(grp_df_prop)[grp_o,]
    rownames(grp) <- NULL
    if (rlang::is_false(remove_ind)) {
      grp$rank <- 1:nrow(grp)
      grp <- grp |> relocate(rank)
    }
    kable(head(grp, n = n_results), format = "simple")
  }

de_seq_run <- function(run_col, dds) {
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$padj),]
  print(resOrdered)
  vsd <- varianceStabilizingTransformation(dds) # vst(dds)
  plotPCA(vsd, intgroup = c(run_col))
}

pre_process_df <- function(df_test) {
  df_test$ERBB4_SEQ <- as.numeric(df_test$ERBB4_SEQ)
  df_test$MDM4_SEQ <- as.numeric(df_test$MDM4_SEQ)
  df_test$LRRN2_SEQ <- as.numeric(df_test$LRRN2_SEQ)
  df_test$PIK3C2B_SEQ <- as.numeric(df_test$PIK3C2B_SEQ)
  df_test[df_test < 0] <- 0
  df_test[is.na(df_test)] <- 0
  df_test <- round(df_test)
  c_ls <- c(
    "Status",
    "ERBB2_SEQ",
    "ERBB2IP_SEQ",
    "ERBB3_SEQ",
    "ERBB4_SEQ",
    "MDM4_SEQ",
    "LRRN2_SEQ",
    "PIK3C2B_SEQ"
  )
  coldata <- data.matrix(df_test[c_ls])
  coldata.T <- t(df_test[c_ls])
  colnames(coldata.T) <- coldata[0,]
  coldata.T <- coldata.T[1:nrow(coldata.T), ]
  cts <- coldata.T
  return(list(countdata = cts, coldata = coldata))
}
