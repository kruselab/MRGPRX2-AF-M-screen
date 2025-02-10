# PACKAGES ---------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(micropan))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(plotROC))
suppressPackageStartupMessages(library(here))


# PDB UTILS --------------------------------------------------------------------

# extract protein sequences from PDB(s); save as multifasta
pdbs_to_fasta = function(dir = "",
                         every_fifth = TRUE, 
                         save_file = TRUE, 
                         prefix = "") {
  pdbs = list.files(path = dir, pattern = "\\.pdb$", full.names = TRUE) 
  if (every_fifth) {
    pdbs = pdbs[seq(1, length(pdbs), 5)]
  }
  output = ""
  count = 1
  cat(paste0("Parsing ", length(pdbs), " files.\n\n"))
  for (pdb_file in pdbs) {
    cat(paste0("Parsing file ", count, ".\n"))
    data = parse_pdb_to_CA(pdb_file)
    chains = extract_chains(data)
    for (i in chains) {
      chain_data = data %>% filter(chain == i)
      seq = str_flatten(names(AMINO_ACID_CODE)[match(chain_data$res_name, toupper(AMINO_ACID_CODE))])
      output = paste0(output, ">", pdb_file, " chain ", i, "\n", seq, "\n\n")
    }
    output = paste0(output, "\n\n")
    count = count + 1
  }
  if (save_file) {
    write_file(output, paste0(prefix, "seqs.fa"))
  }
  return(invisible())
}


# parse PDB and filter on C alphas
parse_pdb_to_CA = function(pdb_file) {
  raw = read_lines(pdb_file)
  calphas = as.data.frame(raw[grepl("^ATOM", raw)]) %>%
    separate_wider_position(col = 1, 
                            widths = c(record_type = 4, 
                                       space1 = 2,
                                       atom_serial = 5, 
                                       space2 = 1, 
                                       atom_name = 4, 
                                       space3 = 1,
                                       res_name = 3, 
                                       space4 = 1,
                                       chain = 1, 
                                       chain_res = 4, 
                                       space5 = 4,
                                       x = 8, 
                                       y = 8, 
                                       z = 8, 
                                       occupancy = 6, 
                                       bfactor = 6, 
                                       segment = 10,
                                       element = 2,
                                       charge = 2)) %>%
    select(-starts_with("space")) %>%
    mutate(across(everything(), str_trim)) %>%
    mutate(across(c("x", "y", "z", "chain_res", "bfactor"), as.numeric)) %>%
    filter(atom_name == "CA") %>%
    rowid_to_column("res")
  return(calphas)
}


# get vector of xyz values from C alpha dataframe and residue index
get_xyz = function(pdb_df, res_i) {
  pdb_df = pdb_df %>% filter(res_i == !!res_i)
  return(c(pdb_df$x, pdb_df$y, pdb_df$z))
}


# extract chain IDs from C alpha dataframe
extract_chains = function(pdb_df) {
  chains = unique(pdb_df$chain)
  return(chains)
}


# calculate centroid of specified C alphas
calculate_centroid = function(pdb_df, chain_res_is, chain_ID) {
  data = pdb_df %>%
    filter(chain == chain_ID) %>%
    filter(chain_res %in% chain_res_is)
  x_centroid = mean(data$x)
  y_centroid = mean(data$y)
  z_centroid = mean(data$z)
  return(c(x_centroid, y_centroid, z_centroid))
}


# calculate distance between two sets of xyz coordinates
dist_xyz = function(coords1, coords2) {
  x1 = coords1[1]
  y1 = coords1[2]
  z1 = coords1[3]
  x2 = coords2[1]
  y2 = coords2[2]
  z2 = coords2[3]
  dist = sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2 + (z1 - z2) ^ 2)
  return(dist)
}


# extract data from pdb (e.g., all-by-all CA distances and pLDDT)
extract_pdb_data = function(pdb_file, 
                            save_file = TRUE, 
                            save_distance_plot = TRUE, 
                            save_plddt_plot = TRUE) {
  pdb_df = parse_pdb_to_CA(pdb_file)
  coords = pdb_df %>% 
    select(res, chain, chain_res, bfactor, x, y, z) %>%
    dplyr::rename(pLDDT = bfactor)
  pdb_df = cross_join(coords, coords, suffix = c("_i", "_j")) %>% 
    right_join(combn(coords$res, 2) %>% 
                 t() %>% 
                 as.data.frame() %>% 
                 setNames(c("res_i", "res_j")), by = join_by(res_i, res_j)) %>%
    rowwise() %>%
    mutate(CA_distance = dist_xyz(c(x_i, y_i, z_i), 
                                  c(x_j, y_j, z_j))) %>%
    ungroup() %>%
    select(-c(x_i, x_j, y_i, y_j, z_i, z_j))
  if (save_file) {
    write_csv(pdb_df, 
              paste0("../", parse_AF_filename(pdb_file), "_PDB_data.csv"))
  }
  if (save_distance_plot) {
    ggsave(paste0(parse_AF_filename(pdb_file), "_CA_distance_plot.pdf"), 
           plot_CA_distance(pdb_df, pdb_file), 
           height = 5, width = 5)
  }
  if (save_plddt_plot) {
    ggsave(paste0(parse_AF_filename(pdb_file), "_pLDDT_plot.pdf"), 
           plot_plddt(pdb_df, pdb_file), 
           height = 2, width = 6)
  }
  return(pdb_df)
}


# plot CA distances
plot_CA_distance = function(pdb_df, pdb_file = "") {
  plot_distance = ggplot(pdb_df, aes(x = res_i, y = res_j, fill = CA_distance)) +
    geom_tile() +
    geom_tile(aes(x = res_j, y = res_i, fill = CA_distance)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), 
                       position = "top", 
                       name = "Residue 1", 
                       breaks = seq(0, 1000, 50)) +
    scale_y_reverse(expand = expansion(mult = c(0, 0)), 
                    name = "Residue 2", 
                    breaks = seq(0, 1000, 50)) +
    scale_fill_distiller(palette = "GnBu", 
                         na.value = brewer.pal(9, "GnBu")[1], 
                         limits = c(0, 25), 
                         breaks = seq(5, 20, 5)) + 
    guides(fill = guide_colorbar(title = expression(paste("C", alpha, "-C", alpha, " (Å)")),
                                 barwidth = 0.5,
                                 barheight = 12,
                                 ticks.colour = "black",
                                 ticks.linewidth = 0.5,
                                 frame.colour = "black",
                                 frame.linewidth = 0.5)) +
    coord_equal() +
    theme_cowplot()
  return(plot_distance)
}


# plot pLDDT
plot_plddt = function(pdb_df, pdb_file = "") {
  plot_plddt = ggplot(pdb_df, aes(x = res_i, y = pLDDT_i, group = chain_i)) +
    geom_line(color = "black", linewidth = 0.25) +
    geom_point(aes(color = pLDDT_i), size = 1, stroke = 0) +
    scale_x_continuous(name = "Residue", 
                       breaks = seq(0, 1000, 50)) +
    scale_y_continuous(name = "pLDDT", 
                       limits = c(0,100), 
                       breaks = seq(0, 100, 20)) +
    scale_color_fermenter(palette = "RdBu", 
                          direction = 1, 
                          breaks = c(0, 50, 70, 90, 100)) +
    guides(color = guide_colorsteps(title = NULL,
                                    barwidth = 0.5,
                                    barheight = 6,
                                    frame.colour = "black",
                                    frame.linewidth = 0.5)) +
    theme_cowplot()
  return(plot_plddt)
}


# PAE UTILS --------------------------------------------------------------------

# parse pAE data to csv
extract_pae_data = function(pae_file, save_file = TRUE, save_plot = TRUE) {
  pae_list = jsonlite::read_json(pae_file, simplifyVector = TRUE)
  pae_df = data.frame(unlist(pae_list[[2]]))
  pae_df = pae_df %>% 
    rowid_to_column("res_i") %>%
    pivot_longer(!res_i, 
                 names_to = "res_j", 
                 names_prefix = "X", 
                 names_transform = as.numeric, 
                 values_to = "pAE")
  pae_df_1 = pae_df %>% 
    filter(res_i < res_j)
  pae_df_2 = pae_df %>% 
    filter(res_i > res_j) %>%
    dplyr::rename(res_i = res_j, res_j = res_i, reciprocal_pAE = pAE)
  pae_df = full_join(pae_df_1, pae_df_2, by = join_by(res_i, res_j)) %>%
    mutate(max_pAE = pae_list[[1]], pTM = pae_list[[3]], ipTM = pae_list[[4]])
  if (save_file) {
    write_csv(pae_df, 
              paste0("../", parse_AF_filename(pae_file), "_pAE_data.csv"))
  }
  if (save_plot) {
    ggsave(paste0(parse_AF_filename(pae_file), "_pAE_plot.pdf"), 
           plot_pae(pae_df, pae_file), 
           height = 5, width = 5)
  }
  return(pae_df)
}


# plot pAE
plot_pae = function(pae_df, pae_file = "") {
  plot_pae = ggplot(pae_df, aes(x = res_j, y = res_i, fill = pAE)) + 
    geom_tile() +
    geom_tile(aes(x = res_i, y = res_j, fill = reciprocal_pAE)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), 
                       position = "top", 
                       name = "Scored residue", 
                       breaks = seq(0, 1000, 50)) +
    scale_y_reverse(expand = expansion(mult = c(0, 0)), 
                    name = "Aligned residue", 
                    breaks = seq(0, 1000, 50)) +
    scale_fill_distiller(palette = "GnBu", 
                         na.value = brewer.pal(9, "GnBu")[1], 
                         limits = c(0, 35), 
                         breaks = seq(5, 30, 5)) + 
    guides(fill = guide_colorbar(title = "PAE (Å)", 
                                 barwidth = 0.5, 
                                 barheight = 12, 
                                 ticks.colour = "black",
                                 ticks.linewidth = 0.5,
                                 frame.colour = "black",
                                 frame.linewidth = 0.5)) +
    coord_equal() +
    theme_cowplot()
  return(plot_pae)
}


# PARSE UTILS ------------------------------------------------------------------

# parse filenames
parse_AF_filename = function(filename) {
  pattern = "^(.*)\\_+\\w*\\_+(.*)\\_+X.*(rank)\\_00(\\d).*(model\\_\\d).*"
  model_info = gsub(pattern, "\\1_\\2_\\3_\\4_\\5", filename)
  return(model_info)
}


# make directory if it doesn't already exist
make_rel_dir = function(parent_dir_path, path_rel_to_parent, new_dir_name) {
  new_path = paste0(parent_dir_path, path_rel_to_parent, new_dir_name)
  if (dir.exists(new_path)) {
    print(paste0("Directory \'", normalizePath(new_path), "\' already exists."))
  } else {
    dir.create(new_path)
    print(paste0("Directory \'", normalizePath(new_path), "\' created."))
  }
  return(normalizePath(new_path))
}


# get full paths for files that match pattern in specified dir
get_file_paths = function(dir, pattern, recurse = FALSE) {
  paths = list.files(path = dir, 
                     pattern = pattern, 
                     full.names = TRUE,
                     recursive = recurse)
  paths = normalizePath(paths)
  return(paths)
}


# compress all files in dir matching pattern
compress_files_in_dir = function(dir, pattern) {
  files = get_file_paths(dir, pattern)
  file_i = 1
  for (file in files) {
    cat(paste0("Compressing file ", 
               file_i, 
               " of ", 
               length(files), 
               ": ", 
               basename(file), 
               "\n"))
    micropan::xzcompress(file, skip = TRUE)
    file_i = file_i + 1
  }
  cat("Compressing complete!\n")
  return(invisible())
}


# uncompress all files in dir matching pattern
uncompress_files_in_dir = function(dir, pattern) {
  files = get_file_paths(dir, pattern)
  file_i = 1
  for (file in files) {
    cat(paste0("Uncompressing file ", 
               file_i, 
               " of ", 
               length(files), 
               ": ", 
               basename(file), 
               "\n"))
    micropan::xzuncompress(file, skip = TRUE)
    file_i = file_i + 1
  }
  cat("Uncompressing complete!\n")
  return(invisible())
}


# remove unnecessary files and folders from AF-M data directory--BE CAREFUL!
clean_data_dir = function(parent_dir) {
  
  # remove unnecessary subdirectories
  folders = list.dirs(path = parent_dir, full.names = TRUE)
  folders = folders[folders != parent_dir]
  num_folders = length(folders)
  unlink(folders, recursive = TRUE)
  cat(paste0(num_folders, " folders deleted.\n"))
  
  # remove unnecessary files
  file_patterns = c("dgram.*\\.xz$",
                    "a3m*\\.xz$",
                    "\\.png$",
                    "\\.txt$",
                    "\\.webp$",
                    "config\\.json",
                    "cite\\.bibtex")
  for (i in 1:length(file_patterns)) {
    files = list.files(path = parent_dir, 
                       pattern = file_patterns[i], 
                       full.names = TRUE)
    num_files = length(folders)
    file.remove(files)
    cat(paste0(num_folders, " files matching pattern \"", file_patterns[i], "\" deleted.\n"))
  }
  
  # count remaining files
  num_files_remaining = length(list.files(path = parent_dir))
  cat(paste0(num_files_remaining, " files remaining.\n"))
  return(invisible())
}


# build models dataframe to link associated AF output files
build_AF_models_df = function(dir, bypass_done_files = FALSE) {
  pdb = data.frame(pdb_file = get_file_paths(dir, "\\.pdb$")) %>%
    mutate(model = parse_AF_filename(basename(pdb_file)))
  pae = data.frame(pae_file = get_file_paths(dir, "scores.*\\.json$")) %>%
    mutate(model = parse_AF_filename(basename(pae_file)))
  models = full_join(pdb, pae, by = join_by(model))
  done_files = get_file_paths(dir, "\\.done\\.txt$")
  if (nrow(models) == 0) {
    stop("No models found in directory.")
  } else if (nrow(models) != 5 * length(done_files)) {
    if (!bypass_done_files) {
      stop("Mismatch between number of model files and .done.txt files.")
    } else {
      warning("Mismatch between number of model files and .done.txt files.",
              immediate. = TRUE)
    }
  }
  return(models)
}


# read pdb or pae data with uncompression and recompression
read_compressed_file = function(file, save_plots = FALSE) {
  uncompressed_file = gsub("(.*)\\.xz$", "\\1", file)
  micropan::xzuncompress(file, skip = TRUE)
  if (str_detect(uncompressed_file, "\\.pdb$")) {
    data = extract_pdb_data(uncompressed_file,
                            save_file = FALSE,
                            save_distance_plot = save_plots,
                            save_plddt_plot = save_plots)
  } else if (str_detect(uncompressed_file, "scores.*\\.json$")) {
    data = extract_pae_data(uncompressed_file,
                            save_file = FALSE,
                            save_plot = save_plots)
  } else {
    micropan::xzcompress(uncompressed_file, skip = TRUE)
    stop("Failed to parse filename")
  }
  micropan::xzcompress(uncompressed_file, skip = TRUE)
  return(data)
}


# read uncompressed pdb or pae data
read_uncompressed_file = function(file, save_plots = FALSE) {
  uncompressed_file = file
  if (str_detect(uncompressed_file, "\\.pdb$")) {
    data = extract_pdb_data(uncompressed_file,
                            save_file = FALSE,
                            save_distance_plot = save_plots,
                            save_plddt_plot = save_plots)
  } else if (str_detect(uncompressed_file, "scores.*\\.json$")) {
    data = extract_pae_data(uncompressed_file,
                            save_file = FALSE,
                            save_plot = save_plots)
  }
  return(data)
}


# parse all AF output files in directory
parse_AF_output_dir = function(dir, save_plots = FALSE, bypass_done_files = FALSE) {
  
  # check for slash at end of dir
  if (!str_detect(dir, "/$")) {
    warning("No forward slash at end of given directory--it has been added.",
            immediate. = TRUE)
    dir = str_c(dir, "/")
  }
  models = build_AF_models_df(dir, bypass_done_files = bypass_done_files)
  output_path = make_rel_dir(dir, "../", "res_pair_data")
  for (i in 1:nrow(models)) {
    print(paste0("Parsing model ", i, " of ", nrow(models), ": ", models$model[i]))
    pdb = read_uncompressed_file(models$pdb_file[i], save_plots = save_plots)
    pae = read_uncompressed_file(models$pae_file[i], save_plots = save_plots)
    full_data = full_join(pdb, pae, by = join_by(res_i, res_j))
    write_csv(full_data, 
              paste0(output_path, "/", models$model[i], "_res_pair_data.csv"))
  }
  print("Parsing complete!")
  return(invisible())
}


# SCORING UTILS ----------------------------------------------------------------

# parse model into three columns in dataframe
add_model_cols = function(data, model_col) {
  data = mutate(data,
                complex = gsub("(.*)\\_rank.*", 
                               "\\1", 
                               basename(!!enquo(model_col))),
                rank = as.integer(gsub(".*rank\\_(\\d).*", 
                                       "\\1", 
                                       !!enquo(model_col))),
                model_num = as.integer(gsub(".*model\\_(\\d).*", 
                                            "\\1", 
                                            !!enquo(model_col))))
  return(data)
}


# read res_pair_data file to dataframe
read_res_pair_file = function(file) {
  data = read_csv(file, show_col_types = FALSE) %>%
    add_model_cols(model_col = file)
  n_chains = length(unique(c(data$chain_i, data$chain_j)))
  if (n_chains != 2) {
    cat(paste0(n_chains, " chains encountered in file: ", basename(file), ". ",
               "Expected number of chains is 2.\n"))
    data = data %>% slice_head(n = 0)
  }
  return(data)
}


# filter contacts (i.e., interchain residue pairs below distance cutoff)
filter_contacts = function(res_pairs_data, distance_cutoff = 2e4, chain = "") {
  contacts = res_pairs_data %>%
    filter(chain_i != chain_j) %>%
    filter(CA_distance <= distance_cutoff)
  if (!str_equal(chain, "")) {
    cat(paste0("Filtering for contacts involving specified chain: ", chain, ".\n"))
    cat(paste0("Unfiltered contacts: ", nrow(contacts), "\n"))
    contacts = contacts %>%
      filter(str_equal(chain, chain_i) | str_equal(chain, chain_j))
    cat(paste0("Filtered contacts: ", nrow(contacts), "\n"))
  }
  return(contacts)
}


# get unique interface residues from contacts_data
get_interface_residues = function(contacts_data) {
  IF_residues = bind_rows(select(contacts_data, 
                                 res = res_i, 
                                 pLDDT = pLDDT_i),
                          select(contacts_data, 
                                 res = res_j, 
                                 pLDDT = pLDDT_j)) %>%
    distinct()
  return(IF_residues)
}


# calculate pDockQ
calc_pdockq = function(avg_IF_plddt, n_contacts) {
  x = avg_IF_plddt * log10(n_contacts)
  pdockq = (0.724 / (1 + exp(-0.052 * (x - 152.611)))) + 0.018
  return(pdockq)
}


# calculate model-level metrics (avg_pAE, avg_pLDDT, etc.)
calc_model_metrics = function(contacts_data) {
  data = contacts_data
  avg_IF_plddt = mean(get_interface_residues(data)$pLDDT)
  metrics = data.frame(complex = data$complex[1],
                       rank = data$rank[1],
                       model_num = data$model_num[1],
                       n_contacts = nrow(data),
                       avg_pAE = mean(c(data$pAE, data$reciprocal_pAE)),
                       avg_pLDDT = avg_IF_plddt,
                       pDockQ = calc_pdockq(avg_IF_plddt, nrow(data)),
                       pTM = mean(data$pTM),
                       ipTM = mean(data$ipTM),
                       rTM = 0.2 * mean(data$pTM) + 0.8 * mean(data$ipTM))
  return(metrics)
}


# calculate model agreement metrics (e.g., avg_model_support)
calc_model_agreement = function(res_pair_data) {
  data = res_pair_data %>%
    group_by(res_i, res_j) %>%
    summarize(n_models = n(), .groups = "keep")
  n_unique_contacts = nrow(data)
  avg_model_support = mean(data$n_models)
  return(c(n_unique_contacts, avg_model_support))
}


# normalize pLDDT to range from 0-1
calc_norm_pLDDT = function(pLDDT) {
  pLDDT_norm = pLDDT / 100
  return(pLDDT_norm)
}


# normalize pAE to range from 0-1
calc_norm_pAE = function(pAE) {
  pAE_norm = (-pAE / 31.75) + 1
  return(pAE_norm)
}


# add combination feature of 6 predictive features to complex_scores df
add_combo_feature_col = function(data) {
  data = mutate(data, combo_feature = best_model_pTM *
                  avg_pTM *
                  calc_norm_pLDDT(best_model_avg_pLDDT) *
                  calc_norm_pLDDT(avg_avg_pLDDT) *
                  calc_norm_pAE(best_model_avg_pAE) *
                  calc_norm_pAE(avg_avg_pAE))
  return(data)
}


# score all res_pair_data.csv files in specified dir
score_res_pairs_dir = function(dir, distance_cutoff = 10) {
  
  # check for slash at end of dir
  if (!str_detect(dir, "/$")) {
    warning("No forward slash at end of given directory--it has been added.",
            immediate. = TRUE)
    dir = str_c(dir, "/")
  }
  
  # build df to associate res_pair files with their complex
  models = data.frame(file = get_file_paths(dir, "res\\_pair\\_data\\.csv")) %>%
    add_model_cols(model_col = file)
  if (nrow(models) == 0) {
    stop("No models found in directory.")
  }
  output_path = make_rel_dir(dir, "../", "score_data")
  
  # iterate through unique complexes and tabulate scores
  model_metrics = data.frame()
  complex_metrics = data.frame()
  unique_complexes = unique(models$complex)
  complex_i = 1
  for (complex in unique_complexes) {
    print(paste0("Scoring complex ", complex_i, " of ", length(unique_complexes), ": ", complex))
    complex_model_metrics = data.frame()
    complex_res_pairs = data.frame()
    complex_models = filter(models, complex == !!complex)
    
    # calculate model-level scores
    for (model in complex_models$file) {
      res_pairs = read_res_pair_file(model) %>%
        filter_contacts(distance_cutoff)
      complex_model_metrics = bind_rows(complex_model_metrics, 
                                        calc_model_metrics(res_pairs))
      complex_res_pairs = bind_rows(complex_res_pairs, 
                                    res_pairs)
    }
    model_metrics = bind_rows(model_metrics, complex_model_metrics)
    
    # calculate complex-level scores
    complex = data.frame(complex = complex_model_metrics$complex[1])
    avg = complex_model_metrics %>%
      select(-c(complex, rank, model_num)) %>%
      summarize(across(everything(), mean)) %>%
      rename_with(~ gsub("(.*)", "avg_\\1", .x))
    best = complex_model_metrics %>%
      filter(rank == 1) %>%
      select(-c(complex, rank, model_num)) %>%
      rename_with(~ gsub("(.*)", "best_model_\\1", .x))
    model_agree = data.frame(n_unique_contacts = calc_model_agreement(complex_res_pairs)[1],
                             avg_model_support = calc_model_agreement(complex_res_pairs)[2])
    row = bind_cols(complex, avg, best, model_agree)
    row = add_combo_feature_col(row)
    complex_metrics = bind_rows(complex_metrics, row)
    complex_i = complex_i + 1
  }
  
  # write outputs
  write_csv(model_metrics, 
            paste0(output_path, "/model_scores.csv"))
  write_csv(complex_metrics, 
            paste0(output_path, "/complex_scores.csv"))
  cat("Scoring complete!\n")
  return(invisible())
}


# MAIN -------------------------------------------------------------------------

run_pipeline = function(dir, 
                        save_model_plots = FALSE, 
                        distance_cutoff = 10,
                        bypass_done_files = FALSE) {
  parse_AF_output_dir(dir = dir, 
                      save_plots = save_model_plots,
                      bypass_done_files = bypass_done_files)
  score_res_pairs_dir(dir = paste0(dir, "../res_pair_data/"), 
                      distance_cutoff = distance_cutoff)
  return(invisible())
}
