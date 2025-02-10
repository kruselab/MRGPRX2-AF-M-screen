# EVALUATION UTILS -------------------------------------------------------------

# combine positive and negative datasets
combine_metrics = function(pos_file, neg_file) {
  pos = read_csv(pos_file, show_col_types = FALSE) %>% mutate(positive = 1)
  neg = read_csv(neg_file, show_col_types = FALSE) %>% mutate(positive = 0)
  merge = bind_rows(pos, neg)
  return(merge)
}


# add positive column to complex score data
add_positive_col = function(file) {
  data = read_csv(file, show_col_types = FALSE) %>% mutate(positive = 1)
  return(data)
}


# reverse signs on specified columns
reverse_col_signs = function(data, cols = c("avg_avg_pAE", "best_model_avg_pAE")) {
  data = mutate(data, across(!!!enquos(cols), ~ .x * -1))
  return(data)
}


# plot all-by-all scatter matrix
plot_all_metric_scatter = function(data,
                                   title = "",
                                   save_plot = FALSE,
                                   prefix = "") {
  
  plot_scatter_matrix = GGally::ggpairs(select(data, -complex), 
                                        aes(color = as.factor(positive)),
                                        diag = list(continuous = wrap("barDiag", 
                                                                      bins = 30, 
                                                                      alpha = 0.5)), 
                                        lower = list(continuous = wrap("points", 
                                                                       size = 1, 
                                                                       stroke = 0, 
                                                                       alpha = 0.5)), 
                                        progress = FALSE) +
    ggtitle(title) +
    theme_cowplot()
  if (save_plot) {
    ggsave(paste0(prefix, "scatter_matrix.pdf"), 
           plot_scatter_matrix, 
           height = 30, width = 30)
  }
  return(plot_scatter_matrix)
}


# plot ROC
plot_ROC = function(data,
                    title = "",
                    save_plot = FALSE,
                    save_AUCs = FALSE,
                    prefix = "") {
  melt_data = melt_roc(data, "positive", 2:18)
  plot_ROC_basic = ggplot(melt_data, aes(d = D.positive, m = M, color = name)) + 
    geom_roc(n.cuts = 0, size = 0.5)
  AUCs = calc_auc(plot_ROC_basic) %>% 
    select(name, AUC)
  plot_ROC = plot_ROC_basic +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = "dotted") +
    scale_x_continuous(name = "FPR") +
    scale_y_continuous(name = "TPR") +
    ggtitle(title) +
    scale_color_discrete(guide = guide_legend("Metric"),
                         labels = paste0(AUCs$name, ": ", round(AUCs$AUC, 4))) +
    coord_equal() +
    theme_cowplot() +
    theme(legend.text.align = 1)
  if (save_plot) {
    ggsave(paste0(prefix, "ROC.pdf"), 
           plot_ROC, 
           height = 6, width = 8)
  }
  if (save_AUCs) {
    write_csv(AUCs, 
              paste0(prefix, "AUCs.csv"))
  }
  return(plot_ROC)
}


# calculate precision from top slice ranked by specified column
calc_precision_from_slice = function(data, col, slice_n) {
  data = slice_max(data, order_by = !!sym(col), n = slice_n, with_ties = FALSE)
  precision = sum(data$positive) / slice_n
  return(precision)
}


# calculate precision for all slice depths on all columns
calc_all_precisions_by_slice = function(data) {
  cols = names(select(data, -c(complex, positive)))
  slices = 1:nrow(data)
  slice_precisions = data.frame()
  for (col in cols) {
    for (slice in slices) {
      row = data.frame(metric = col,
                       n = slice,
                       precision = calc_precision_from_slice(data, col, slice))
      slice_precisions = bind_rows(slice_precisions, row)
    }
  }
  return(slice_precisions)
}


# plot precision by test (aka slice) depth
plot_precision_by_slice = function(data,
                             title = "",
                             save_plot = FALSE,
                             save_precisions = FALSE,
                             prefix = "") {
  
  # get basic data
  n_models = nrow(data)
  data = calc_all_precisions_by_slice(data)
  if (save_precisions) {
    write_csv(data,
              paste0(prefix, "precisions_by_slice.csv"))
  }
  
  # rank order metrics by precision at 5% depth
  rank_order = data %>%
    group_by(metric) %>%
    filter(n <= n_models / 20) %>%
    summarize(mean_precision_top_slice = mean(precision)) %>%
    arrange(desc(mean_precision_top_slice)) %>%
    mutate(metric_factor = factor(metric, levels = metric))
  data = mutate(data, metric_factor = factor(metric, levels = rank_order$metric))
  
  # plot
  plot_precision = ggplot(data, aes(x = n, y = precision)) +
    geom_rect(xmin = 0, 
              xmax = n_models / 20, 
              ymin = 0, 
              ymax = 1, 
              fill = "gray80", 
              alpha = 0.1) +
    geom_line() +
    geom_hline(data = rank_order,
               aes(yintercept = mean_precision_top_slice),
               color = "gray80", 
               linetype = "dotted") + 
    geom_text(data = rank_order, 
              aes(label = round(mean_precision_top_slice, 3)),
              hjust = 1,
              vjust = 1,
              x = n_models,
              y = 1) +
    scale_x_continuous(name = "Test depth (%)",
                       breaks = seq(0, n_models, length.out = 5),
                       labels = seq(0, 100, length.out = 5)) +
    scale_y_continuous(name = "True Positive Rate", 
                       expand = expansion(mult = c(0, 0.05))) +
    ggtitle(title) +
    facet_wrap(vars(metric_factor)) +
    coord_fixed(n_models) +
    theme_cowplot()
  if (save_plot) {
    ggsave(paste0(prefix, "precision_by_slice.pdf"), 
           plot_precision, 
           height = 10, width = 10)
  } 
  return(plot_precision)
}


# run evaluation pipeline on pos and neg complex_metrics data
evaluate_complex_metrics = function(pos_file, 
                                    neg_file, 
                                    plot_title = "", 
                                    prefix = "") {
  data = combine_metrics(pos_file, neg_file) %>%
    reverse_col_signs()
  plot_ROC(data,
           title = plot_title, 
           save_plot = TRUE, 
           save_AUCs = TRUE,
           prefix = prefix)
  plot_all_metric_scatter(data, 
                          title = plot_title, 
                          save_plot = TRUE, 
                          prefix = prefix)
  plot_precision_by_slice(data,
                    title = plot_title, 
                    save_plot = TRUE, 
                    save_precisions = TRUE,
                    prefix = prefix)
  return(invisible())
}


# CLONE SELECTION UTILS --------------------------------------------------------

# combine similar tables using recursive file search
combine_tables = function(dir, pattern, prefix = "") {
  files = get_file_paths(dir = dir, 
                         pattern = pattern,
                         recurse = TRUE)
  output = data.frame()
  i = 1
  for (file in files) {
    cat(paste0("Processing file ", i, " of ", length(files), ": ", basename(file), "\n"))
    data = read_csv(file, show_col_types = FALSE)
    output = bind_rows(output, data)
    i = i + 1
  }
  write_csv(output, 
            paste0(prefix, "combined_table.csv"))
  return(invisible())
}


# plot histograms of a "discovery mode" dataset alongside TP and TN datasets
plot_discovery_validation_compare = function(discovery_file,
                                             discovery_name,
                                             TP_file,
                                             TN_file,
                                             save_plot = TRUE) {
  # read and merge data
  TP = read_csv(TP_file, show_col_types = FALSE) %>% 
    mutate(dataset = "TP")
  TN = read_csv(TN_file, show_col_types = FALSE) %>% 
    mutate(dataset = "TN")
  discovery = read_csv(discovery_file, show_col_types = FALSE) %>% 
    mutate(dataset = discovery_name)
  merge = bind_rows(TP, TN, discovery)
  
  # hit rate calculations
  n_above_TN = sum(discovery$combo_feature > max(TN$combo_feature))
  n_total = nrow(discovery)
  frac_above_TN = round(100 * n_above_TN / n_total, 3)
  hit_rate_message = paste0(
    n_above_TN,
    "/",
    n_total,
    " (",
    frac_above_TN,
    "%) clones have higher score than maximum in TN set.")
  cat(hit_rate_message)
  
  # plot
  plot = ggplot(data = merge, aes(x = combo_feature, 
                                  y = after_stat(density),
                                  group = dataset)) +
    labs(title = paste0(discovery_name, " vs validation set"),
         subtitle = hit_rate_message) +
    geom_histogram(binwidth = 0.01, 
                   boundary = 0, 
                   fill = "#1f78b4") +
    geom_vline(xintercept = max(TN$combo_feature), 
               linetype = "dashed", 
               color = "gray80") +
    scale_x_continuous("Combination feature") +
    scale_y_continuous("Density") +
    facet_wrap(facets = vars(dataset), ncol = 1) +
    theme_cowplot()
  
  # save plot
  disc_name_no_whitespace = gsub("\\s", "\\_", discovery_name)
  if (save_plot) {
    ggsave(paste0(disc_name_no_whitespace, "_validation_comparison.pdf"), 
           plot, 
           height = 6, 
           width = 6)
  }
  return(plot)
}

