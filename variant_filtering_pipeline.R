## =========================
## Variant filtering pipeline
## =========================


## Homozygous filtering

## 0) Libraries
library(dplyr)
library(stringr)
library(readr)
library(openxlsx)
library(readxl)

## =========================
## 0) Input tables for parent lookup
## =========================

lgc_path <- "C:\\Users\\Sherry\\OneDrive - University College London\\Desktop\\Internship in Queen Square\\Rauan_LGC.xlsx"
# Set the path to Rauan_LGC file

# Read proband parent RK table from Rauan_LGC folder

map_path <- "C:\\Users\\Sherry\\OneDrive - University College London\\Desktop\\Internship in Queen Square\\combined_sequencing_file_locations.xlsx"
# Read RK to LI mapping table

stopifnot(file.exists(lgc_path))
stopifnot(file.exists(map_path))

## Helper: normalise RK to RK0001 format
norm_rk <- function(x) {
  x <- toupper(as.character(x))
  rk <- stringr::str_extract(x, "RK\\d+")
  n <- suppressWarnings(as.integer(stringr::str_remove(rk, "^RK")))
  ifelse(is.na(n), NA_character_, sprintf("RK%04d", n))
}

##extract RK ID
rk_num <- function(rk) {
  suppressWarnings(as.integer(stringr::str_remove(rk, "^RK")))
}


## -------------------------
## A) Read RK to LI mapping
## -------------------------
map_raw <- readxl::read_xlsx(map_path, sheet = 1)

## IMPORTANT FIX:
## In this project, LI numeric id is stored in LogID and needs "LI" prefix
## RK is extracted from File_location or File_data_id
rk_to_li <- map_raw %>%   
  dplyr::mutate(
    li_id = ifelse(
      is.na(LogID) | as.character(LogID) == "",
      NA_character_,
      paste0("LI", as.character(as.integer(LogID)))
    ),
    # LI id derived from LogID column
    
    rk_raw = dplyr::coalesce(
      stringr::str_extract(as.character(File_location), "RK\\d+"),
      stringr::str_extract(as.character(File_data_id), "RK\\d+")
    ),
    # Extract RK from mapping file from file location and file data id columns
    
    rk = norm_rk(rk_raw)
    # Normalise RK format
  ) %>%
  dplyr::filter(!is.na(li_id), li_id != "", !is.na(rk), rk != "") %>%
  dplyr::distinct(rk, li_id)
## From the table produced by mutate(), retain only rows where both li_id and rk are present and using distinct to remove duplicate rows, keeping only one

## Also create LI to RK mapping (folder first strategy)
li_to_rk <- rk_to_li %>% dplyr::distinct(li_id, rk)
# LI to RK lookup table

## -------------------------
## B)Read LGC and build proband to parent RK mapping using proband based row grouping
## -------------------------
wb_lgc <- openxlsx::loadWorkbook(lgc_path)
lgc_sheet <- "Main sheet"
# Use Main sheet explicitly because it contains Participate.type, Clinical.details

lgc_raw <- openxlsx::readWorkbook(wb_lgc, sheet = lgc_sheet)
# Read LGC main sheet

## Column detection
dna_name  <- grep("^DNA(\\.|\\s)number$", names(lgc_raw), ignore.case = TRUE, value = TRUE)[1]
clin_name <- grep("^Clinical(\\.|\\s)details$", names(lgc_raw), ignore.case = TRUE, value = TRUE)[1]
part_name <- grep("^Participate(\\.|\\s)type$", names(lgc_raw), ignore.case = TRUE, value = TRUE)[1]


stopifnot(!is.na(dna_name))
stopifnot(!is.na(clin_name))
stopifnot(!is.na(part_name))

## Keyword patterns
father_pat  <- "\\b(father|dad|paternal)\\b"
mother_pat  <- "\\b(mother|mum|maternal)\\b"
proband_pat <- "proband"
# Role detection patterns (Participate.type preferred, Clinical.details fallback)

lgc_keep <- lgc_raw %>%
  dplyr::mutate(
    dna_number = as.character(.data[[dna_name]]),
    clinical_details = as.character(.data[[clin_name]]),
    participate_type = as.character(.data[[part_name]]),
    
    rk_raw = stringr::str_extract(dna_number, "RK\\d+"),
    rk = norm_rk(rk_raw),
    
    pt_low = stringr::str_to_lower(stringr::str_trim(participate_type)),
    ##remove whitespace then convert to lowercase
    
    rel = dplyr::case_when(
      stringr::str_detect(pt_low, father_pat) ~ "father",
      stringr::str_detect(pt_low, mother_pat) ~ "mother",
      stringr::str_detect(pt_low, proband_pat) ~ "proband",
      TRUE ~ NA_character_
    ),
    
    is_proband = !is.na(rel) & rel == "proband",
    row_index = dplyr::row_number(),
    proband_group = cumsum(is_proband),
    proband_group = dplyr::if_else(proband_group == 0, NA_integer_, proband_group)
  ) %>%
  dplyr::filter(
    !is.na(rk), rk != "",
    !is.na(proband_group)
  )
# Keep only rows with valid RK and assigned proband_group

cat("\nLGC role counts (raw):\n")
print(table(lgc_keep$rel, useNA = "ifany"))
# Sanity check

## Build pedigree from proband based grouping (order dependent on row sequence)
pedig_rk <- lgc_keep %>%
  dplyr::group_by(proband_group) %>%
  dplyr::summarise(
    proband_row = row_index[which(rel == "proband")[1]],
    proband_rk = rk[which(rel == "proband")[1]],
    father_rk  = {
      idx <- which(rel == "father" & row_index > proband_row)
      if (length(idx) == 0) NA_character_ else rk[idx[1]]
    },
    mother_rk  = {
      idx <- which(rel == "mother" & row_index > proband_row)
      if (length(idx) == 0) NA_character_ else rk[idx[1]]
    },
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(proband_rk), proband_rk != "")
# One inferred trio per proband_group

cat("\nPedigree rows built (proband groups with proband): ", nrow(pedig_rk), "\n")
# Debug: how many families are usable

## -------------------------
## C) Folder first master builder (LI in folder -> RK -> proband_group pedigree -> parents -> LI)
## -------------------------
build_master_for_run <- function(folder_li_vec) {
  
  folder_li_vec <- unique(na.omit(folder_li_vec))
  # LI list actually present in this folder
  
  ## 1) Map folder LI to RK
  folder_li_rk <- li_to_rk %>%
    dplyr::filter(li_id %in% folder_li_vec) %>%
    dplyr::distinct(li_id, rk)
  
  if (nrow(folder_li_rk) == 0) {
    return(tibble::tibble(
      proband_group = integer(),
      proband_rk = character(),
      father_rk = character(),
      mother_rk = character(),
      proband_li = character(),
      father_li = character(),
      mother_li = character(),
      status = character()
    ))
  }
  
  ## 2) Use the RK range represented in this folder
  rk_min <- min(rk_num(folder_li_rk$rk), na.rm = TRUE)
  rk_max <- max(rk_num(folder_li_rk$rk), na.rm = TRUE)
  
  ## 3) Keep only inferred probands within this RK range
  proband_master <- pedig_rk %>%
    dplyr::mutate(proband_rk_num = rk_num(proband_rk)) %>%
    dplyr::filter(!is.na(proband_rk_num), proband_rk_num >= rk_min, proband_rk_num <= rk_max) %>%
    dplyr::transmute(
      proband_group = proband_group,
      proband_rk = proband_rk,
      father_rk = father_rk,
      mother_rk = mother_rk
    ) %>%
    dplyr::distinct()
  
  ## 4) Map proband RK to LI (may be missing)
  proband_master <- proband_master %>%
    dplyr::left_join(rk_to_li, by = c("proband_rk" = "rk")) %>%
    dplyr::rename(proband_li = li_id)
  
  ## 5) Map parent RK to LI
  proband_master <- proband_master %>%
    dplyr::left_join(rk_to_li, by = c("father_rk" = "rk")) %>%
    dplyr::rename(father_li = li_id) %>%
    dplyr::left_join(rk_to_li, by = c("mother_rk" = "rk")) %>%
    dplyr::rename(mother_li = li_id)
  
  ## 6) Add run status
  proband_master <- proband_master %>%
    dplyr::mutate(
      status = dplyr::case_when(
        is.na(proband_li) | proband_li == "" ~ "no LI ID",
        !(proband_li %in% folder_li_vec) ~ "no LI file",
        TRUE ~ "ok"
      )
    ) %>%
    dplyr::mutate(proband_rk_num = rk_num(proband_rk)) %>%
    dplyr::arrange(proband_rk_num) %>%
    dplyr::select(-proband_rk_num)
  
  proband_master
}
## =========================
## 1) Show remaining rows + first n rows
## =========================
show_state <- function(df, step_name, n = 10) {
  cat("\n====================\n")
  cat(step_name, "\n")
  cat("Remaining rows:", nrow(df), "\n")
  cat("====================\n")
  
  print(
    df %>%
      select(chrom, pos, gene, gt, impact, consequence, count_het, count_hom, ad, polyphen,
             gnomad_af, gnomadg_af_afr, cadd_phred) %>%
      head(n),
    width = Inf
  )
}


## =========================
## 2) Helpers for phenotype and gene grouping
## =========================

## Keep only consecutive duplicated rows in a given column
keep_consecutive_duplicates <- function(x) {
  x <- str_trim(as.character(x))
  x == dplyr::lag(x) | x == dplyr::lead(x)
}

## Create a run_id for consecutive blocks of identical values
make_run_id <- function(x) {
  x2 <- str_trim(as.character(x))
  cumsum(x2 != dplyr::lag(x2, default = x2[1])) + 1L
}

## Show state for phenotype or symbol branches (auto handle group columns)
show_state_pheno_symbol <- function(df, step_name, n = 10) {
  cat("\n====================\n")
  cat(step_name, "\n")
  cat("Remaining rows:", nrow(df), "\n")
  cat("====================\n")
  
  base_cols <- c("chrom", "pos", "gene", "gt", "impact", "symbol", "consequence",
                 "count_het", "count_hom", "ad",
                 "gnomad_af", "gnomadg_af_afr", "cadd_phred", "omim_pheno")
  
  group_cols <- c("group_id", "group_source", "group_value", "run_id")
  
  cols_to_show <- c(base_cols, intersect(group_cols, names(df)))
  
  print(
    df %>%
      select(all_of(cols_to_show)) %>%
      head(n)
  )
}

## Keep only groups where all rows pass a condition
group_all_pass <- function(df, pass_col) {
  df %>%
    group_by(group_id) %>%
    filter(all({{ pass_col }})) %>%
    ungroup()
}

## =========================
## 3) Clean and standardise types (so filtering works consistently)
## =========================
clean_df <- function(df) {
  
  ## Force commonly inconsistent annotation columns to character
  to_char_cols <- c("protein_position", "cdna_position", "cds_position",
                    "amino_acids", "codons")
  
  for (cc in intersect(to_char_cols, names(df))) {
    df[[cc]] <- as.character(df[[cc]])
  }
  
  df %>%
    mutate(
      ## -------------------------
      ## Standardise chrom and pos
      ## -------------------------
      chrom = str_trim(as.character(chrom)),
      
      ## Extract main chromosome number/letter
      chrom_main = str_to_upper(chrom),
      chrom_main = str_replace(chrom_main, "^CHR", ""),
      chrom_main = str_extract(chrom_main, "^(\\d+|X|Y|M|MT)"),
      chrom_main = dplyr::case_when(
        chrom_main %in% c("M") ~ "MT",
        TRUE ~ chrom_main
      ),
      # Make chrom consistent: 1-22, X, Y, MT
      
      chrom = chrom_main,
      # IMPORTANT FIX: use cleaned chromosome for parent check
      
      pos = suppressWarnings(as.integer(readr::parse_number(as.character(pos)))),
      
      # Make pos consistent integer
      
      ## -------------------------
      ## Keep original cleaning
      ## -------------------------
      gt = as.character(gt),
      impact = as.character(impact),
      consequence = as.character(consequence),
      
      count_het = str_trim(as.character(count_het)),
      count_hom = str_trim(as.character(count_hom)),
      
      ad = str_trim(as.character(ad)),
      sift = str_trim(as.character(sift)),
      polyphen = str_trim(as.character(polyphen)),
      
      gnomad_af = readr::parse_number(as.character(gnomad_af)),
      gnomadg_af_afr = readr::parse_number(as.character(gnomadg_af_afr)),
      cadd_phred = readr::parse_number(as.character(cadd_phred))
    ) %>%
    mutate(
      ad_clean = str_remove_all(ad, "\\[|\\]"),
      ad_ref = suppressWarnings(as.integer(str_extract(ad_clean, "^\\s*\\d+"))),
      ad_alt = suppressWarnings(as.integer(str_extract(ad_clean, "(?<=,)\\s*\\d+")))
    )
}


## =========================
## 3A) Parent comparison helpers
## =========================

GT_HET <- c("0_1", "1_0", "0|1", "1|0")
# Heterozygous genotype encodings

get_parent_gt_at_pos <- function(parent_df, chrom_in, pos_in) {
  # Return parent gt at chrom and pos, NA if not found or no parent df
  
  if (is.null(parent_df) || nrow(parent_df) == 0) return(NA_character_)
  
  chrom_in <- as.character(chrom_in)
  pos_in <- suppressWarnings(as.integer(pos_in))
  
  hit <- parent_df %>%
    dplyr::filter(.data$chrom == chrom_in, .data$pos == pos_in) %>%
    dplyr::slice(1)
  
  if (nrow(hit) == 0) return(NA_character_)
  
  as.character(hit$gt[1])
}

pos_present_in_parent <- function(parent_df, chrom_in, pos_in) {
  # Return TRUE if parent has any record at chrom and pos
  !is.na(get_parent_gt_at_pos(parent_df, chrom_in, pos_in))
}



## =========================
## 4) Single sample runner (original pipeline, wrapped)
## =========================
run_one_sample <- function(df,
                           sample_id = "sample",
                           father_df = NULL,
                           mother_df = NULL,
                           verbose = FALSE,
                           proband_li = NA_character_,
                           proband_rk = NA_character_) {
  
  # Run filtering for one sample and optionally check parents
  
  ## Define genotype groups
  GT_HOM_ALT <- c("1_1", "1|1")
  GT_HET     <- c("0_1", "1_0", "0|1", "1|0")
  
  ## Quick check
  stopifnot(is.data.frame(df) || tibble::is_tibble(df))
  
  ## Clean / standardise
  df0 <- clean_df(df)
  orig_cols <- names(df)
  
  if (verbose) {
    show_state(df0, paste0(sample_id, " START: after cleaning"))
  }
  
  ## =========================
  ## Homozygous filtering steps
  ## =========================
  
  ## Step 1: GT keep 1_1 and 1|1 only
  d1 <- df0 %>% filter(gt %in% GT_HOM_ALT)
  if (verbose) show_state(d1, paste0(sample_id, " HOM Step 1: gt (homozygous)"))
  
  ## Step 2: IMPACT remove MODIFIER and NA, keep others
  d2 <- d1 %>% filter(!is.na(impact), impact != "MODIFIER")
  if (verbose) show_state(d2, paste0(sample_id, " Step 2: impact (remove MODIFIER + NA)"))
  
  ## Step 3: CONSEQUENCE remove synonymous_variant and synonymous_variant&NMD_transcript_variant
  ## (Do NOT remove NA / blank)
  d3 <- d2 %>% filter(
    is.na(consequence) |
      consequence == "" |
      !consequence %in% c("synonymous_variant", "synonymous_variant&NMD_transcript_variant")
  )
  if (verbose) show_state(d3, paste0(sample_id, " Step 3: consequence (remove synonymous types, keep NA/blank)"))
  
  ## Step 4: count_het keep 0-10 and blank/NA
  d4 <- d3 %>% filter(
    is.na(count_het) | count_het == "" | count_het %in% as.character(0:10)
  )
  if (verbose) show_state(d4, paste0(sample_id, " Step 4: count_het (keep 0-10 + blank/NA)"))
  
  ## Step 5: count_hom keep 0-5 and blank/NA
  d5 <- d4 %>% filter(
    is.na(count_hom) | count_hom == "" | count_hom %in% as.character(0:5)
  )
  if (verbose) show_state(d5, paste0(sample_id, " Step 5: count_hom (keep 0-5 + blank/NA)"))
  
  ## Step 6: AD remove rows where (first number > 0) OR (second number < 10)
  ## Meaning: keep only when ad_ref <= 0 AND ad_alt >= 10
  d6 <- d5 %>% filter(!(ad_ref > 0 | ad_alt < 10))
  if (verbose) show_state(d6, paste0(sample_id, " Step 6: AD (remove ad_ref>0 OR ad_alt<10)"))
  
  ## Step 7: remove rows where sift or polyphen contains "benign"
  d7 <- d6 %>% filter(
    (is.na(sift) | !str_detect(tolower(sift), "benign")) &
      (is.na(polyphen) | !str_detect(tolower(polyphen), "benign"))
  )
  if (verbose) show_state(d7, paste0(sample_id, " Step 7: sift + polyphen (remove benign)"))
  
  ## Step 8: gnomad_af remove > 0.001
  d8 <- d7 %>% filter(
    is.na(gnomad_af) | gnomad_af <= 0.001
  )
  if (verbose) show_state(d8, paste0(sample_id, " Step 8: gnomad_af (<= 0.001)"))
  
  ## Step 9: gnomadg_af_afr remove > 0.05
  d9 <- d8 %>% filter(
    is.na(gnomadg_af_afr) | gnomadg_af_afr <= 0.05
  )
  if (verbose) show_state(d9, paste0(sample_id, " Step 9: gnomadg_af_afr (<= 0.05)"))
  
  ## Step 10: cadd_phred remove < 15
  d10 <- d9 %>% filter(
    is.na(cadd_phred) | cadd_phred >= 15
  )
  if (verbose) show_state(d10, paste0(sample_id, " Step 10: cadd_phred (>= 15) FINAL"))
  
  ## Summary table of counts per step (HOM)
  counts_hom <- tibble(
    sample = sample_id,
    stream = "Homozygous",
    step = c("df0", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9", "d10"),
    n = c(nrow(df0), nrow(d1), nrow(d2), nrow(d3), nrow(d4), nrow(d5),
          nrow(d6), nrow(d7), nrow(d8), nrow(d9), nrow(d10))
  )
  
  ## =========================
  ## Heterozygous filtering pipeline
  ## ========================= 
  
  ## Step 1: GT keep 1|0, 0|1, 1_0, 0_1
  h1 <- df0 %>% filter(gt %in% GT_HET)
  if (verbose) show_state(h1, paste0(sample_id, " HET Step 1: gt (heterozygous)"))
  
  ## Step 2: IMPACT remove MODIFIER and NA
  h2 <- h1 %>% filter(!is.na(impact), impact != "MODIFIER")
  if (verbose) show_state(h2, paste0(sample_id, " Step 2: impact (remove MODIFIER+NA)"))
  
  ## Step 3: CONSEQUENCE remove synonymous_variant and synonymous_variant&NMD_transcript_variant
  ## (Do NOT remove NA / blank)
  h3 <- h2 %>% filter(
    is.na(consequence) |
      consequence == "" |
      !consequence %in% c("synonymous_variant",
                          "synonymous_variant&NMD_transcript_variant")
  )
  if (verbose) show_state(h3, paste0(sample_id, " Step 3: consequence (remove synonymous types)"))
  
  ## Step 4: count_het keep 0-10 and blank/NA
  h4 <- h3 %>% filter(
    is.na(count_het) | count_het == "" | count_het %in% as.character(0:10)
  )
  if (verbose) show_state(h4, paste0(sample_id, " Step 4: count_het (keep 0-10 + blank/NA)"))
  
  ## Step 5: count_hom keep 0-5 and blank/NA
  h5 <- h4 %>% filter(
    is.na(count_hom) | count_hom == "" | count_hom %in% as.character(0:5)
  )
  if (verbose) show_state(h5, paste0(sample_id, " Step 5: count_hom (keep 0-5 + blank/NA)"))
  
  ## Step 6: remove rows where sift or polyphen contains "benign"
  h6 <- h5 %>% filter(
    (is.na(sift) | !str_detect(tolower(sift), "benign")) &
      (is.na(polyphen) | !str_detect(tolower(polyphen), "benign"))
  )
  if (verbose) show_state(h6, paste0(sample_id, " HET Step 6: sift + polyphen (remove benign)"))
  
  ## =========================
  ## Step 6: phenotype / gene based filtering (two branch strategy)
  ## =========================
  
  ## -------------------------
  ## Branch A:
  ## Keep variants with consecutive duplicated phenotype descriptions
  ## in omim_pheno column (must be adjacent and identical)
  ## -------------------------
  
  A_pheno <- h6 %>%
    mutate(omim_pheno = str_trim(as.character(omim_pheno))) %>%   # clean phenotype text
    filter(
      !is.na(omim_pheno),                                         # remove NA
      omim_pheno != "",                                           # remove blank
      keep_consecutive_duplicates(pick(omim_pheno)$omim_pheno)     # keep only consecutive identical phenotypes
    ) %>%
    mutate(
      group_source = "A_pheno",                                   # group comes from phenotype branch
      group_value = omim_pheno,                                   # group defined by phenotype text
      run_id = make_run_id(group_value),                          # consecutive block id
      group_id = paste(group_source, group_value, run_id, sep = " | ")
    )
  
  show_state_pheno_symbol(A_pheno, "Branch A: consecutive duplicated omim_pheno")
  
  ## -------------------------
  ## Branch B step 1:
  ## Keep variants with blank / NA omim_pheno
  ## -------------------------
  
  B_blank <- h6 %>%
    mutate(omim_pheno = str_trim(as.character(omim_pheno))) %>%
    filter(is.na(omim_pheno) | omim_pheno == "")
  
  show_state_pheno_symbol(B_blank, "Branch B-1: omim_pheno blank only")
  
  ## -------------------------
  ## Branch B step 2:
  ## Within blank phenotype variants, keep only consecutive duplicated genes (symbol)
  ## -------------------------
  
  B_gene <- B_blank %>%
    mutate(symbol = str_trim(as.character(symbol))) %>%
    filter(
      !is.na(symbol),
      symbol != "",
      keep_consecutive_duplicates(pick(symbol)$symbol)
    ) %>%
    mutate(
      group_source = "B_symbol",
      group_value = symbol,
      run_id = make_run_id(group_value),
      group_id = paste(group_source, group_value, run_id, sep = " | ")
    )
  
  show_state_pheno_symbol(B_gene, "Branch B-2: consecutive duplicated symbol within blank phenotype")
  
  ## Step 6 final:
  ## Merge Branch A and Branch B, keep groups together by sorting
  h7 <- bind_rows(A_pheno, B_gene) %>%
    arrange(group_source, group_value, run_id)
  if (verbose) show_state_pheno_symbol(h7, paste0(sample_id, " Step 6: merged phenotype + gene prioritisation (grouped)"))
  
  ## Step 7: AD filter for heterozygous (group based)
  ## Rule: ad_ref >= 10, ad_alt >= 10, and abs(ad_ref - ad_alt) <= 30
  h8 <- h7 %>%
    mutate(
      pass_ad = !is.na(ad_ref) & !is.na(ad_alt) &
        ad_ref >= 10 &
        ad_alt >= 10 &
        abs(ad_ref - ad_alt) <= 30
    ) %>%
    group_all_pass(pass_ad)
  if (verbose) show_state_pheno_symbol(h8, paste0(sample_id, " HET Step 7: AD (both >=10, diff <=30)"))
  
  ## Step 8: gnomad_af remove > 0.001 (group based)
  h9 <- h8 %>%
    mutate(pass_gnomad_af = is.na(gnomad_af) | gnomad_af <= 0.001) %>%
    group_all_pass(pass_gnomad_af)
  if (verbose) show_state_pheno_symbol(h9, paste0(sample_id, " HET Step 8: gnomad_af (<= 0.001)"))
  
  ## Step 9: gnomadg_af_afr remove > 0.05 (group based)
  h10 <- h9 %>%
    mutate(pass_afr = is.na(gnomadg_af_afr) | gnomadg_af_afr <= 0.05) %>%
    group_all_pass(pass_afr)
  if (verbose) show_state_pheno_symbol(h10, paste0(sample_id, " HET Step 9: gnomadg_af_afr (<= 0.05)"))
  
  ## Step 10: cadd_phred remove < 15 (group based)
  h11 <- h10 %>%
    mutate(pass_cadd = is.na(cadd_phred) | cadd_phred >= 15) %>%
    group_all_pass(pass_cadd)
  if (verbose) show_state_pheno_symbol(h11, paste0(sample_id, " HET Step 10: cadd_phred (>= 15) FINAL"))
  
  ## =========================
  ## Parent check after filtering
  ## =========================
  
  hom_reason <- NA_character_
  het_reason <- NA_character_
  # Store parent check status separately for HOM and HET
  
  ## Clean parent dfs using the same cleaning function
  father0 <- if (!is.null(father_df)) clean_df(father_df) else NULL
  # Prepare father data for position search
  
  mother0 <- if (!is.null(mother_df)) clean_df(mother_df) else NULL
  # Prepare mother data for position search
  
  
  ## -------------------------
  ## HOM parent check
  ## Rule: require both parents and both parents must be HET at same position
  ## -------------------------
  if (nrow(d10) > 0) {
    
    if (is.null(father0) && is.null(mother0)) {
      d10 <- d10[0, ]
      hom_reason <- "no parent data"
    } else if (is.null(father0)) {
      d10 <- d10[0, ]
      hom_reason <- "no father file"
    } else if (is.null(mother0)) {
      d10 <- d10[0, ]
      hom_reason <- "no mother file"
    } else {
      
      d10_before <- nrow(d10)
      
      d10 <- d10 %>%
        rowwise() %>%
        mutate(
          f_gt = get_parent_gt_at_pos(father0, chrom, pos),
          m_gt = get_parent_gt_at_pos(mother0, chrom, pos),
          pass_parents = (!is.na(f_gt) && f_gt %in% GT_HET) &&
            (!is.na(m_gt) && m_gt %in% GT_HET)
        ) %>%
        ungroup() %>%
        filter(pass_parents) %>%
        select(-f_gt, -m_gt, -pass_parents)
      
      if (nrow(d10) == 0) {
        hom_reason <- "failed parent check"
      } else {
        hom_reason <- "ok"
      }
    }
    
  } else {
    hom_reason <- "no variants after filtering"
  }
  # If no HOM variants, still report a reason
  
  ## =========================
  ## HET parent check ON (compound heterozygous segregation)
  ## Rule:
  ##   For each group_id, require exactly 2 variants.
  ##   Keep group if variants segregate across parents:
  ##     one variant present in father only and the other present in mother only
  ##   If only one parent is available, keep group only if exactly one of the 2 variants is present in that parent.
  ## =========================
  
  if (nrow(h11) > 0) {
    
    ## Clean parent dfs already prepared as father0 / mother0 above
    
    if (is.null(father0) && is.null(mother0)) {
      h11 <- h11[0, ]
      het_reason <- "no parent data"
    } else {
      
      ## Keep only groups with exactly 2 variants (compound het must be 2 hits)
      h11 <- h11 %>%
        group_by(group_id) %>%
        filter(n() == 2) %>%
        ungroup()
      
      if (nrow(h11) == 0) {
        het_reason <- "failed parent check (not exactly 2 variants per group)"
      } else {
        
        ## Decide which groups pass segregation
        h11_keep <- h11 %>%
          group_by(group_id) %>%
          summarise(
            chrom1 = chrom[1],
            pos1   = pos[1],
            chrom2 = chrom[2],
            pos2   = pos[2],
            
            ## Presence means the exact chrom+pos exists in parent file
            f1 = if (!is.null(father0)) pos_present_in_parent(father0, chrom1, pos1) else NA,
            f2 = if (!is.null(father0)) pos_present_in_parent(father0, chrom2, pos2) else NA,
            m1 = if (!is.null(mother0)) pos_present_in_parent(mother0, chrom1, pos1) else NA,
            m2 = if (!is.null(mother0)) pos_present_in_parent(mother0, chrom2, pos2) else NA,
            
            ## Keep if segregating across parents (or weak check if single parent)
            keep_group = {
              if (!is.null(father0) && !is.null(mother0)) {
                (f1 && !m1 && m2 && !f2) || (f2 && !m2 && m1 && !f1)
              } else if (!is.null(father0) && is.null(mother0)) {
                (f1 != f2)
              } else if (is.null(father0) && !is.null(mother0)) {
                (m1 != m2)
              } else {
                FALSE
              }
            },
            .groups = "drop"
          ) %>%
          filter(keep_group) %>%
          select(group_id)
        
        ## Keep only passing groups
        h11 <- h11 %>% semi_join(h11_keep, by = "group_id")
        
        if (nrow(h11) == 0) {
          het_reason <- "failed parent check"
        } else {
          het_reason <- "ok"
        }
      }
    }
    
  } else {
    het_reason <- "no variants after filtering"
  }
  
  
  
  ## -------------------------
  ## Combine into one display reason
  ## -------------------------
  parent_reason <- paste0("HOM: ", hom_reason, " | HET: ", het_reason)
  # Used for Excel id row
  
  ## Summary table of counts per step (HET)
  counts_het <- tibble(
    sample = sample_id,
    stream = "Heterozygous",
    step = c("df0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", "h9", "h10","h11"),
    n = c(nrow(df0), nrow(h1), nrow(h2), nrow(h3), nrow(h4), nrow(h5),
          nrow(h6), nrow(h7), nrow(h8), nrow(h9), nrow(h10),nrow(h11))
  )
  
  ## =========================
  ## Export results to Excel (FULL COLUMNS, HOM + blank + HET)
  ## =========================
  
  ## 1) Prepare HOM output (keep all columns)
  hom_out <- d10 %>% mutate(stream = "Homozygous")
  
  ## 2) Prepare HET output (keep grouped order, keep all columns)
  het_out <- h11 %>%
    arrange(group_source, group_value, run_id) %>%
    mutate(stream = "Heterozygous")
  
  ## 3) Align columns safely (do not error if some columns do not exist)
  common_cols <- union(names(hom_out), names(het_out))
  hom_out2 <- hom_out %>% select(any_of(common_cols))
  het_out2 <- het_out %>% select(any_of(common_cols))
  
  ## 4) Create one blank row in the middle (same columns as combined output)
  blank_row <- tibble::as_tibble(
    as.list(setNames(rep(list(NA), length(common_cols)), common_cols))
  )
  
  ## 5) Combine: HOM -> blank row -> HET
  combined <- bind_rows(hom_out2, blank_row, het_out2)
  
  ## =========================
  ## Build display table for Excel (human readable)
  ## =========================
  
  ## Helper: insert a blank row after each HET group
  insert_blank_between_groups <- function(df, group_col, blank_row) {
    if (nrow(df) == 0) return(df)
    
    out <- list()
    groups <- unique(df[[group_col]])
    
    for (g in groups) {
      out[[length(out) + 1]] <- df %>% dplyr::filter(.data[[group_col]] == g)
      out[[length(out) + 1]] <- blank_row
    }
    
    dplyr::bind_rows(out)
  }
  
  ## Clean proband ID (remove trailing _annotated if present)
  clean_id <- sub("annotated$", "", sample_id)
  clean_id <- sub("_$", "", clean_id)
  
  ## Template based on original patient columns only
  base_template <- df0 %>%
    dplyr::select(dplyr::any_of(orig_cols)) %>%
    dplyr::slice(0)
  
  ## Add display only columns (stream is used for yellow highlighting later)
  base_template <- base_template %>%
    dplyr::mutate(
      proband_id = NA_character_,
      stream = NA_character_,
      group_id = NA_character_
    ) %>%
    dplyr::relocate(proband_id, stream, group_id, .before = 1)
  
  ## HOM rows
  hom_disp <- hom_out2 %>%
    dplyr::select(dplyr::any_of(orig_cols)) %>%
    dplyr::mutate(
      proband_id = NA_character_,
      stream = "HOM",
      group_id = NA_character_
    ) %>%
    dplyr::relocate(proband_id, stream, group_id, .before = 1)
  
  ## HET rows (keep group_id temporarily for grouping)
  het_disp <- het_out2 %>%
    dplyr::select(dplyr::any_of(c(orig_cols, "group_id"))) %>%
    dplyr::mutate(
      proband_id = NA_character_,
      stream = "HET"
    ) %>%
    dplyr::relocate(proband_id, stream, group_id, .before = 1)
  
  ## One blank row (type safe)
  blank_row_display <- base_template
  blank_row_display[1, ] <- NA
  
  ## Proband ID row
  id_row <- blank_row_display
  # Create proband id row
  
  ## Prefer LI from master, fallback to clean_id
  show_li <- if (!is.na(proband_li) && proband_li != "") proband_li else clean_id
  
  ## Attach RK if available
  show_id <- if (!is.na(proband_rk) && proband_rk != "") {
    paste0(show_li, " | ", proband_rk)
  } else {
    show_li
  }
  
  if (!is.na(parent_reason) && parent_reason != "") {
    id_row$proband_id <- paste0(show_id, " | ", parent_reason)
  } else {
    id_row$proband_id <- show_id
  }
  
  # Append negative reason next to proband id if needed
  
  ## Negative row
  neg_row <- blank_row_display
  neg_row$proband_id <- "negative"
  
  ## Insert blank rows between HET groups
  het_disp_grouped <- het_disp %>%
    dplyr::arrange(group_id) %>%
    insert_blank_between_groups("group_id", blank_row_display)
  
  ## Does this proband have any results
  has_results <- nrow(hom_out2) > 0 || nrow(het_out2) > 0
  
  ## Final display
  if (!has_results) {
    display_table <- dplyr::bind_rows(
      id_row,
      blank_row_display,
      neg_row,
      blank_row_display,
      blank_row_display
    )
  } else {
    display_table <- dplyr::bind_rows(
      id_row,
      blank_row_display,
      hom_disp,
      blank_row_display,
      het_disp_grouped,
      blank_row_display,
      blank_row_display
    )
  }
  
  ## Keep only proband_id + stream + original patient columns (drop group_id)
  display_table <- display_table %>%
    dplyr::select(proband_id, stream, dplyr::any_of(orig_cols))
  
  cat("\n[DEBUG] sample:", sample_id, "\n")
  cat("[DEBUG] HET reason:", het_reason, " | HET rows after parent check:", nrow(h11), "\n")
  
  
  ## 6) Return outputs
  list(
    sample = sample_id,
    Results = combined,
    Homozygous = hom_out2,
    Heterozygous = het_out2,
    Display = display_table,
    Counts = bind_rows(counts_hom, counts_het)
  )
}


## =========================
## 5) Batch IO
## =========================

## Read one annotated file into df
read_annotated_file <- function(path, sheet = 1) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("xlsx", "xls", "xlsm")) {
    readxl::read_excel(path, sheet = sheet) %>% as.data.frame()
  } else if (ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE) %>% as.data.frame()
  } else if (ext == "tsv") {
    readr::read_tsv(path, show_col_types = FALSE) %>% as.data.frame()
  } else if (ext == "rds") {
    readRDS(path)
  } else {
    stop("Unsupported file type: ", ext)
  }
}

## Make valid Excel sheet names (max 31 chars, remove invalid characters)
safe_sheet_name <- function(x) {
  x2 <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  x2 <- gsub("\\s+", " ", x2)
  x2 <- trimws(x2)
  substr(x2, 1, 31)
}

make_note_block <- function(proband_li = NA_character_,
                            proband_rk = NA_character_,
                            note = "",
                            orig_cols = character()) {
  
  if (length(orig_cols) == 0) {
    base <- tibble::tibble(proband_id = NA_character_, stream = NA_character_)
  } else {
    base <- tibble::as_tibble(
      as.list(setNames(rep(list(NA), length(orig_cols)), orig_cols))
    ) %>%
      dplyr::mutate(
        proband_id = NA_character_,
        stream = NA_character_
      ) %>%
      dplyr::relocate(proband_id, stream)
  }
  
  id_row <- base
  label <- if (!is.na(proband_li) && proband_li != "") {
    paste0(proband_li, " | ", proband_rk, " | ", note)
  } else {
    paste0(proband_rk, " | ", note)
  }
  id_row$proband_id <- label
  
  neg_row <- base
  neg_row$proband_id <- "negative"
  
  blank_row <- base
  blank_row[1, ] <- NA
  
  dplyr::bind_rows(
    id_row,
    blank_row,
    neg_row,
    blank_row,
    blank_row
  )
}

## Run multiple files and write all outputs into one Excel workbook
run_batch_to_one_excel <- function(in_dir,
                                   out_xlsx,
                                   pattern = "\\.(xlsx|xls|csv|tsv|rds)$",
                                   excel_sheet = 1,
                                   verbose = FALSE) {
  
  # List all input files
  paths <- list.files(in_dir, full.names = TRUE, pattern = "\\.(csv|tsv|xlsx|xlsm|xls|rds)$", ignore.case = TRUE)
  
  # Stop if no files
  if (length(paths) == 0) {
    stop("No input files found in: ", in_dir)
  }
  
  # Extract LI IDs from filenames in this folder
  sample_ids_all <- tools::file_path_sans_ext(basename(paths))
  li_in_files_all <- stringr::str_extract(sample_ids_all, "LI\\d+")
  
  # Build master for all probands, then annotate whether LI/file exists in this folder
  master_folder <- build_master_for_run(li_in_files_all)
  
  if (nrow(master_folder) == 0) {
    stop("No probands found via LGC.")
  }
  
  cat("Run master built. Probands found:", nrow(master_folder), "\n")
  print(table(master_folder$status, useNA = "ifany"))
  
  # Prepare template original columns for note blocks
  template_orig_cols <- character()
  first_existing_file <- paths[1]
  if (length(first_existing_file) > 0 && !is.na(first_existing_file)) {
    tmp_df <- read_annotated_file(first_existing_file, sheet = excel_sheet)
    template_orig_cols <- names(tmp_df)
  }
  
  # Create output workbook
  wb <- openxlsx::createWorkbook()
  all_display <- list()
  
  # Run probands in RK order
  run_list <- master_folder %>%
    dplyr::mutate(proband_rk_num = rk_num(proband_rk)) %>%
    dplyr::arrange(proband_rk_num) %>%
    dplyr::select(-proband_rk_num)
  
  # Loop over each proband row
  for (i in seq_len(nrow(run_list))) {
    
    parent_row <- run_list[i, , drop = FALSE]
    
    proband_li <- parent_row$proband_li
    proband_rk <- parent_row$proband_rk
    status <- parent_row$status
    
    # Case 1: no LI ID
    if (status == "no LI ID") {
      all_display[[paste0("RK_", proband_rk)]] <- make_note_block(
        proband_li = NA_character_,
        proband_rk = proband_rk,
        note = "no LI ID",
        orig_cols = template_orig_cols
      )
      next
    }
    
    # Case 2: no LI file
    if (status == "no LI file") {
      all_display[[paste0("RK_", proband_rk)]] <- make_note_block(
        proband_li = proband_li,
        proband_rk = proband_rk,
        note = "no LI file",
        orig_cols = template_orig_cols
      )
      next
    }
    
    # Case 3: file exists, run full analysis
    p <- paths[stringr::str_detect(basename(paths), stringr::fixed(proband_li))][1]
    
    if (length(p) == 0 || is.na(p)) {
      all_display[[paste0("RK_", proband_rk)]] <- make_note_block(
        proband_li = proband_li,
        proband_rk = proband_rk,
        note = "no LI file",
        orig_cols = template_orig_cols
      )
      next
    }
    
    # Read proband file
    sample_id <- tools::file_path_sans_ext(basename(p))
    df <- read_annotated_file(p, sheet = excel_sheet)
    
    # Lookup parent LI IDs
    father_li <- parent_row$father_li
    mother_li <- parent_row$mother_li
    
    # Find father file path
    father_path <- if (!is.na(father_li) && father_li != "") {
      hit <- paths[stringr::str_detect(basename(paths), stringr::fixed(father_li))]
      if (length(hit) > 0) hit[1] else NA_character_
    } else NA_character_
    
    # Find mother file path
    mother_path <- if (!is.na(mother_li) && mother_li != "") {
      hit <- paths[stringr::str_detect(basename(paths), stringr::fixed(mother_li))]
      if (length(hit) > 0) hit[1] else NA_character_
    } else NA_character_
    
    # Read parent files (if exist)
    father_df <- if (!is.na(father_path)) read_annotated_file(father_path, sheet = excel_sheet) else NULL
    mother_df <- if (!is.na(mother_path)) read_annotated_file(mother_path, sheet = excel_sheet) else NULL
    
    # Run full filtering + parent check
    res <- run_one_sample(
      df,
      sample_id = sample_id,
      father_df = father_df,
      mother_df = mother_df,
      verbose = verbose,
      proband_li = proband_li,
      proband_rk = proband_rk
    )
    
    # Store display table
    all_display[[paste0("RK_", proband_rk)]] <- res$Display
  }
  
  # Stop if no probands processed
  if (length(all_display) == 0) {
    stop("No probands available for output.")
  }
  
  # Combine all proband display outputs
  final_table <- dplyr::bind_rows(all_display)
  
  ## =========================
  ## GLOBAL sort proband blocks by RK number
  ## =========================
  
  extract_rk_from_proband_id <- function(x) {
    rk <- stringr::str_extract(as.character(x), "RK\\d+")
    rk_num(rk)
  }
  
  id_idx <- which(!is.na(final_table$proband_id) & stringr::str_detect(final_table$proband_id, "RK\\d+"))
  
  if (length(id_idx) > 0) {
    start_idx <- id_idx
    end_idx <- c(id_idx[-1] - 1, nrow(final_table))
    
    blocks_tbl <- lapply(seq_along(start_idx), function(i) {
      final_table[start_idx[i]:end_idx[i], , drop = FALSE]
    })
    
    blocks_rk <- vapply(blocks_tbl, function(b) {
      extract_rk_from_proband_id(b$proband_id[1])
    }, numeric(1))
    
    ord <- order(is.na(blocks_rk), blocks_rk)
    final_table <- dplyr::bind_rows(blocks_tbl[ord])
  }
  
  # Write final table
  openxlsx::addWorksheet(wb, "All_probands")
  openxlsx::writeData(wb, "All_probands", final_table)
  
  # Highlight HET rows
  het_style <- openxlsx::createStyle(fgFill = "#FFF59D")
  het_rows <- which(final_table$stream == "HET")
  
  if (length(het_rows) > 0) {
    openxlsx::addStyle(
      wb,
      sheet = "All_probands",
      style = het_style,
      rows = het_rows + 1,
      cols = 1:ncol(final_table),
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  # Hide stream column
  stream_col <- which(names(final_table) == "stream")
  if (length(stream_col) == 1) {
    openxlsx::setColWidths(wb, "All_probands", cols = stream_col, widths = 0)
  }
  
  # Save workbook
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  cat("Saved: ", out_xlsx, "\n")
  invisible(out_xlsx)
}

## =========================
## 6) Run batch
## =========================
in_dir <- "C:/Users/Sherry/OneDrive - University College London/Desktop/Internship in Queen Square/Georgia LI ID"

out_xlsx <- file.path(
  "C:/Users/Sherry/OneDrive - University College London/Desktop/Internship in Queen Square",
  paste0("Variant_filtering_ALL_probands_parentcheck_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
)

run_batch_to_one_excel(
  in_dir = in_dir,
  out_xlsx = out_xlsx,
  verbose = FALSE
)
