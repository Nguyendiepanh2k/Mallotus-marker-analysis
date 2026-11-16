library(ape)
library(tidyverse)
library(writexl)

# Set your working directory here
# Fill in your FASTA file paths below (one alignment per marker, aligned sequences)

seqs <- list(
  rbcL  = read.FASTA(),
  matK  = read.FASTA(),
  ITS   = read.FASTA(),
  phyC  = read.FASTA(),
  ncpGS = read.FASTA(),
  ITS_phyC      = read.FASTA(),
  matK_phyC     = read.FASTA(),
  rbcL_matK_ITS = read.FASTA()
)
lapply(seqs, length)

# Adjust species_lookup to your own reference names.

species_lookup <- c(
  "Mape_ref" = "Mape",
  "Mbar_ref" = "Mbar",
  "Mpan_ref" = "Mpan",
  "Mpel_ref" = "Mpel",
  "Mphi_ref" = "Mphi",
  "Mflo_ref" = "Mflo"
)

extract_species <- function(seq_names, lookup) {
  sapply(seq_names, function(name) {
    if (grepl("^[A-Z][a-z]+_[a-z]+", name)) {
      if (name %in% names(lookup)) return(lookup[name])
      else return(name)
    } else {
      return(substr(name, 1, 4))
    }
  })
}

species_list <- lapply(seqs, function(s) extract_species(names(s), species_lookup))
lapply(species_list, table)

dist_list <- lapply(seqs, function(s) {
  dist.dna(s, model = "raw", pairwise.deletion = TRUE)
})

marker_names <- names(dist_list)

get_per_species_vals <- function(dist_matrix, species) {  
  mat     <- as.matrix(dist_matrix)  
  seq_nms <- rownames(mat)  
  sp_list <- unique(species[species != "OtherSp"])
  
  setNames(lapply(sp_list, function(sp) {  
    focal_idx <- which(species == sp)  
    other_idx <- which(species != sp)
    
    intra_vals <- c()  
    if (length(focal_idx) > 1)  
      for (i in 1:(length(focal_idx) - 1))  
        for (j in (i + 1):length(focal_idx)) {  
          val <- mat[focal_idx[i], focal_idx[j]]  
          if (!is.na(val)) intra_vals <- c(intra_vals, val)  
        }
    
    inter_vals <- c(); inter_keys <- character(0)  
    for (i in focal_idx)  
      for (j in other_idx) {  
        val <- mat[i, j]  
        if (!is.na(val)) {  
          inter_vals <- c(inter_vals, val)  
          inter_keys <- c(inter_keys,  
                          paste(sort(c(seq_nms[i], seq_nms[j])), collapse = "__"))  
        }  
      }  
    list(intra = intra_vals, inter = setNames(inter_vals, inter_keys))  
  }), sp_list)  
}  

per_species_vals <- setNames(
  lapply(marker_names, function(m)
    get_per_species_vals(dist_list[[m]], species_list[[m]])),
  marker_names
)

intra_list <- lapply(per_species_vals, function(pv) unlist(lapply(pv, `[[`, "intra")))
inter_list <- lapply(per_species_vals, function(pv) {  
  v <- unlist(lapply(pv, `[[`, "inter"))  
  v[!duplicated(names(v))]         
}) 
stopifnot(all(sapply(inter_list, function(v) !any(duplicated(names(v))))))
lapply(intra_list, length)
lapply(inter_list, length)

summary_stats <- do.call(rbind, lapply(marker_names, function(marker) {
  intra <- intra_list[[marker]]
  inter <- inter_list[[marker]]
  data.frame(
    marker      = marker,
    mean_intra  = round(mean(intra, na.rm = TRUE), 4),
    sd_intra    = round(sd(intra,   na.rm = TRUE), 4),
    min_intra   = round(min(intra,  na.rm = TRUE), 4),
    max_intra   = round(max(intra,  na.rm = TRUE), 4),
    mean_inter  = round(mean(inter, na.rm = TRUE), 4),
    sd_inter    = round(sd(inter,   na.rm = TRUE), 4),
    min_inter   = round(min(inter,  na.rm = TRUE), 4),
    max_inter   = round(max(inter,  na.rm = TRUE), 4),
    barcode_gap = round(mean(inter, na.rm = TRUE) - mean(intra, na.rm = TRUE), 4)
  )
}))
print(summary_stats)

wilcox_intra_inter <- do.call(rbind, lapply(marker_names, function(marker) {
  intra <- intra_list[[marker]][!is.na(intra_list[[marker]])]
  inter <- inter_list[[marker]][!is.na(inter_list[[marker]])]
  n_intra <- length(intra)
  n_inter <- length(inter)
  test <- wilcox.test(inter, intra, paired = FALSE, alternative = "greater")
  w_plus  <- test$statistic
  w_minus <- (as.numeric(n_inter) * as.numeric(n_intra)) - w_plus
  
  data.frame(
    marker      = marker,
    W_plus      = w_plus,
    W_minus     = w_minus,
    p_value     = round(test$p.value, 5),
    significant = ifelse(test$p.value < 0.05, "Yes", "No")
  )
}))
print(wilcox_intra_inter)

# Per species × marker:
# overlap_count = number of inter values < max(intra)
# (inter values that "fall into" the intra range → bad for discrimination)
# Marker-level overlap_pct (POOLED):
# = sum(overlap_count across all species) / sum(total_inter across all species) × 100
# This is the per-species incidence, NOT the average of per-species percentages.

overlap_by_species <- do.call(rbind, lapply(marker_names, function(marker) {
  pv <- per_species_vals[[marker]]
  do.call(rbind, lapply(names(pv), function(sp) {
    intra_v <- pv[[sp]]$intra
    inter_v <- pv[[sp]]$inter
    
    if (length(intra_v) == 0 || length(inter_v) == 0) {
      return(data.frame(
        marker = marker, species = sp,
        n_intra = length(intra_v), n_inter = length(inter_v),
        max_intra = NA, min_inter = NA,
        overlap_exists = NA, overlap_count = NA,
        total_inter = NA, overlap_pct = NA
      ))
    }
    
    max_intra   <- max(intra_v, na.rm = TRUE)
    min_inter   <- min(inter_v, na.rm = TRUE)
    ovlp_cnt    <- sum(inter_v < max_intra, na.rm = TRUE)
    total_inter <- sum(!is.na(inter_v))
    
    data.frame(
      marker         = marker,
      species        = sp,
      n_intra        = length(intra_v),
      n_inter        = length(inter_v),
      max_intra      = round(max_intra, 4),
      min_inter      = round(min_inter, 4),
      overlap_exists = ifelse(min_inter < max_intra, "Yes", "No"),
      overlap_count  = ovlp_cnt,
      total_inter    = total_inter,
      overlap_pct    = round(ovlp_cnt / total_inter * 100, 2)
    )
  }))
}))

overlap_summary <- overlap_by_species %>%
  filter(!is.na(overlap_count)) %>%
  group_by(marker) %>%
  summarise(
    n_species_evaluated      = n(),
    n_species_with_overlap   = sum(overlap_exists == "Yes"),
    pct_species_with_overlap = round(n_species_with_overlap / n_species_evaluated * 100, 1),
    total_overlap_count      = sum(overlap_count),       # pooled numerator
    total_inter_count        = sum(total_inter),         # pooled denominator
    overlap_pct_pooled       = round(total_overlap_count / total_inter_count * 100, 2),
    .groups = "drop"
  )

print(overlap_by_species)
print(overlap_summary)

# Mnan is absent from phyC and ncpGS, so it is excluded from the pairwise cross-marker Wilcoxon tests.

get_named_dists <- function(dist_matrix, species, exclude_sp = "Mnan") {
  mat     <- as.matrix(dist_matrix)
  seq_nms <- rownames(mat)
  
  keep    <- which(species != exclude_sp)
  mat     <- mat[keep, keep]
  seq_nms <- seq_nms[keep]
  species <- species[keep]
  
  sp_list <- unique(species[species != "OtherSp"])
  
  intra_named <- numeric(0)
  inter_named <- numeric(0)
  
  for (sp in sp_list) {
    focal_idx <- which(species == sp)
    other_idx <- which(species != sp)
    
    if (length(focal_idx) > 1) {
      for (i in 1:(length(focal_idx) - 1))
        for (j in (i + 1):length(focal_idx)) {
          val <- mat[focal_idx[i], focal_idx[j]]
          if (!is.na(val)) {
            key <- paste(sort(c(seq_nms[focal_idx[i]], seq_nms[focal_idx[j]])),
                         collapse = "__")
            intra_named[key] <- val
          }
        }
    }
    
    for (i in focal_idx)
      for (j in other_idx) {
        val <- mat[i, j]
        if (!is.na(val)) {
          key <- paste(sort(c(seq_nms[i], seq_nms[j])), collapse = "__")
          inter_named[key] <- val
        }
      }
  }
  
  list(intra = intra_named, inter = inter_named)
}

named_dists <- setNames(
  lapply(marker_names, function(m)
    get_named_dists(dist_list[[m]], species_list[[m]], exclude_sp = "Mnan")),
  marker_names
)

cat("\n── Common intra pairs between markers ──\n")
for (i in 1:(length(marker_names) - 1))
  for (j in (i + 1):length(marker_names))
    cat(sprintf("  %s vs %s : %d pairs\n", marker_names[i], marker_names[j],
                length(intersect(names(named_dists[[marker_names[i]]]$intra),
                                 names(named_dists[[marker_names[j]]]$intra)))))

cat("\n── Common inter pairs between markers ──\n")
for (i in 1:(length(marker_names) - 1))
  for (j in (i + 1):length(marker_names))
    cat(sprintf("  %s vs %s : %d pairs\n", marker_names[i], marker_names[j],
                length(intersect(names(named_dists[[marker_names[i]]]$inter),
                                 names(named_dists[[marker_names[j]]]$inter)))))

run_pairwise_wilcox_pairs <- function(dist_type, alternative) {
  do.call(rbind, lapply(1:(length(marker_names) - 1), function(i) {
    do.call(rbind, lapply((i + 1):length(marker_names), function(j) {
      m1 <- marker_names[i]; m2 <- marker_names[j]
      v1 <- named_dists[[m1]][[dist_type]]
      v2 <- named_dists[[m2]][[dist_type]]
      
      common_keys <- intersect(names(v1), names(v2))
      x <- v1[common_keys]; y <- v2[common_keys]
      valid <- !is.na(x) & !is.na(y)
      x <- x[valid]; y <- y[valid]
      
      if (length(x) < 3) {
        return(data.frame(
          marker_1 = m1, marker_2 = m2, n_pairs = length(x),
          mean_m1 = NA, mean_m2 = NA, V = NA,
          p_m1_wins = NA, p_m2_wins = NA, winner = "Insufficient data"
        ))
      }
      
      alt_rev   <- ifelse(alternative == "less", "greater", "less")
      t_m1_wins <- wilcox.test(x, y, paired = TRUE, alternative = alternative)
      t_m2_wins <- wilcox.test(x, y, paired = TRUE, alternative = alt_rev)
      
      data.frame(
        marker_1  = m1,
        marker_2  = m2,
        n_pairs   = length(x),
        mean_m1   = round(mean(x), 4),
        mean_m2   = round(mean(y), 4),
        V         = t_m1_wins$statistic,
        p_m1_wins = round(t_m1_wins$p.value, 5),
        p_m2_wins = round(t_m2_wins$p.value, 5),
        winner    = ifelse(t_m1_wins$p.value < 0.05, m1,
                           ifelse(t_m2_wins$p.value < 0.05, m2, "No difference"))
      )
    }))
  }))
}

pairwise_intra <- run_pairwise_wilcox_pairs("intra", "less")
names(pairwise_intra)[names(pairwise_intra) == "mean_m1"] <- "mean_intra_m1"
names(pairwise_intra)[names(pairwise_intra) == "mean_m2"] <- "mean_intra_m2"
names(pairwise_intra)[names(pairwise_intra) == "p_m1_wins"] <- "p_m1_lower_intra"
names(pairwise_intra)[names(pairwise_intra) == "p_m2_wins"] <- "p_m2_lower_intra"
names(pairwise_intra)[names(pairwise_intra) == "winner"]    <- "winner_lower_intra"
print(pairwise_intra)

pairwise_inter <- run_pairwise_wilcox_pairs("inter", "greater")
names(pairwise_inter)[names(pairwise_inter) == "mean_m1"] <- "mean_inter_m1"
names(pairwise_inter)[names(pairwise_inter) == "mean_m2"] <- "mean_inter_m2"
names(pairwise_inter)[names(pairwise_inter) == "p_m1_wins"] <- "p_m1_higher_inter"
names(pairwise_inter)[names(pairwise_inter) == "p_m2_wins"] <- "p_m2_higher_inter"
names(pairwise_inter)[names(pairwise_inter) == "winner"]    <- "winner_higher_inter"
print(pairwise_inter)

per_species_summary <- setNames(
  lapply(marker_names, function(marker) {
    pv <- per_species_vals[[marker]]
    do.call(rbind, lapply(names(pv), function(sp) {
      intra_v    <- pv[[sp]]$intra
      inter_v    <- pv[[sp]]$inter
      mean_intra <- ifelse(length(intra_v) > 0, mean(intra_v, na.rm = TRUE), NA)
      mean_inter <- mean(inter_v, na.rm = TRUE)
      data.frame(
        species     = sp,
        mean_intra  = round(mean_intra, 4),
        mean_inter  = round(mean_inter, 4),
        barcode_gap = round(mean_inter - ifelse(is.na(mean_intra), 0, mean_intra), 4)
      )
    }))
  }),
  marker_names
)

per_species_common <- lapply(per_species_summary, function(df) {
  df <- df[df$species != "Mnan", ]
  df[order(df$species), ]
})
lapply(per_species_common, nrow)
lapply(per_species_common, function(df) df$species)

mrg <- merge(per_species_common[[m1]][, c("species","barcode_gap")],  
             per_species_common[[m2]][, c("species","barcode_gap")],  
             by = "species", suffixes = c(".1", ".2"))  
x <- mrg$barcode_gap.1; y <- mrg$barcode_gap.2  
valid <- !is.na(x) & !is.na(y); x <- x[valid]; y <- y[valid]  

pairwise_gap <- do.call(rbind, lapply(1:(length(marker_names) - 1), function(i) {
  do.call(rbind, lapply((i + 1):length(marker_names), function(j) {
    m1 <- marker_names[i]; m2 <- marker_names[j]
    x  <- per_species_common[[m1]]$barcode_gap
    y  <- per_species_common[[m2]]$barcode_gap
    valid <- !is.na(x) & !is.na(y)
    x <- x[valid]; y <- y[valid]
    
    t_m1_wins <- wilcox.test(x, y, paired = TRUE, alternative = "greater")
    t_m2_wins <- wilcox.test(x, y, paired = TRUE, alternative = "less")
    
    data.frame(
      marker_1          = m1,
      marker_2          = m2,
      n_species         = length(x),
      mean_gap_m1       = round(mean(x), 4),
      mean_gap_m2       = round(mean(y), 4),
      V                 = t_m1_wins$statistic,
      p_m1_larger_gap   = round(t_m1_wins$p.value, 5),
      p_m2_larger_gap   = round(t_m2_wins$p.value, 5),
      winner_larger_gap = ifelse(t_m1_wins$p.value < 0.05, m1,
                                 ifelse(t_m2_wins$p.value < 0.05, m2, "No difference"))
    )
  }))
}))
print(pairwise_gap)

# Win = p < 0.05 (one-sided Wilcoxon) in the correct direction.
# wins_intra = times this marker had significantly SMALLER intra
# wins_inter = times this marker had significantly LARGER  inter
# wins_gap   = times this marker had significantly LARGER  gap
# composite_score = wins_intra + wins_inter + wins_gap

count_wins <- function(pw_df, winner_col) {
  all_m <- unique(c(pw_df$marker_1, pw_df$marker_2))
  setNames(
    sapply(all_m, function(m) sum(pw_df[[winner_col]] == m, na.rm = TRUE)),
    all_m
  )
}

intra_wins <- count_wins(pairwise_intra, "winner_lower_intra")
inter_wins <- count_wins(pairwise_inter, "winner_higher_inter")
gap_wins   <- count_wins(pairwise_gap,   "winner_larger_gap")

ranking <- summary_stats %>%
  left_join(
    overlap_summary %>%
      select(marker, total_overlap_count, total_inter_count, overlap_pct_pooled,
             n_species_with_overlap, pct_species_with_overlap),
    by = "marker"
  ) %>%
  left_join(
    wilcox_intra_inter %>%
      select(marker, W_plus, W_minus, p_ranksum = p_value, sig_ranksum = significant),
    by = "marker"
  ) %>%
  mutate(
    wins_intra      = as.integer(intra_wins[marker]),
    wins_inter      = as.integer(inter_wins[marker]),
    wins_gap        = as.integer(gap_wins[marker]),
    composite_score = wins_intra + wins_inter + wins_gap
  ) %>%
  arrange(
    n_species_with_overlap,    
    total_overlap_count,    
    desc(barcode_gap),        
    mean_intra,               
    desc(mean_inter)          
  ) %>%
  mutate(rank = row_number()) %>%
  select(rank, marker,
         n_species_with_overlap, total_overlap_count, overlap_pct_pooled,
         barcode_gap, mean_intra, sd_intra, mean_inter, sd_inter,
         W_plus, W_minus, p_ranksum, sig_ranksum,
         composite_score)                 

print(ranking)

df_all <- do.call(rbind, lapply(marker_names, function(marker) {
  data.frame(
    marker   = marker,
    type     = c(rep("Intraspecific", length(intra_list[[marker]])),
                 rep("Interspecific", length(inter_list[[marker]]))),
    distance = c(intra_list[[marker]], inter_list[[marker]])
  )
}))

write_xlsx(
list(
  "Summary Statistics"      = summary_stats,
  "Rank-Sum Inter vs Intra" = wilcox_intra_inter,
  "Overlap by Species"      = overlap_by_species,
  "Overlap Summary"         = as.data.frame(overlap_summary),
  "Pairwise Intra"          = pairwise_intra,
  "Pairwise Inter"          = pairwise_inter,
  "Pairwise Barcode Gap"    = pairwise_gap,
  "Ranking"                 = ranking,
  "All Distances"           = df_all
),
path = "barcode_results.xlsx"
)

plot_boxplot <- ggplot(df_all, aes(x = marker, y = distance, fill = type)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = c("Intraspecific" = "#4C72B0",
                               "Interspecific" = "#DD8452")) +
  labs(title = "Intra- vs Interspecific P-Distances by Marker",
       x = "Marker", y = "P-Distance", fill = "") +
  theme_classic(base_size = 13)

plot_gap <- ggplot(ranking, aes(x = reorder(marker, barcode_gap),
                                y = barcode_gap, fill = barcode_gap)) +
  geom_col(show.legend = FALSE) +
  scale_fill_gradient(low = "#f7c59f", high = "#c0392b") +
  labs(title = "Barcode Gap by Marker",
       x = "Marker", y = "Interspecific − Intraspecific (P-Distance)") +
  coord_flip() +
  theme_classic(base_size = 13)

plot_overlap <- ggplot(
  overlap_by_species %>% filter(!is.na(overlap_pct)),
  aes(x = marker, y = overlap_pct, fill = marker)
) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, show.legend = FALSE) +
  geom_point(data = as.data.frame(overlap_summary),
             aes(x = marker, y = overlap_pct_pooled),
             shape = 18, size = 3, color = "black") +
  labs(title = "Overlap % per Marker",
       subtitle = "Boxplot = per-species distribution  |  ◆ = pooled overlap %",
       x = "Marker", y = "Overlap % (inter values < max intra)") +
  theme_classic(base_size = 13)

print(plot_boxplot)
print(plot_gap)
print(plot_overlap)

ggsave("boxplot_distances.svg", plot = plot_boxplot, width = 8, height = 6)
ggsave("barcode_gap.svg",       plot = plot_gap,     width = 8, height = 6)
ggsave("overlap_pct.svg",       plot = plot_overlap, width = 8, height = 6)
