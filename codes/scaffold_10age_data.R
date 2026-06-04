################################################################################
# SCAFFOLDING SCRIPT — produces 10 age band data files from the original
# Testing 10 age band model
################################################################################

input_dir  <- file.path("data")
output_dir <- file.path("data")

old_ages <- c("0 to 4","5 to 11","12 to 17","18 to 29","30 to 39",
              "40 to 49","50 to 59","60 to 69", "70+")
new_ages <- c("0 to 4","5 to 11","12 to 17","18 to 29","30 to 39",
              "40 to 49","50 to 59","60 to 64","65 to 74","75+")

## --- 1. Demographics --------------------------------------------------------
demog <- read.csv(file.path(input_dir, "demographics2021.csv"),
                  header = TRUE, stringsAsFactors = FALSE)

# For each (IMD, rural) block, expand 9 age rows into 10
blocks <- split(demog, list(demog$IMD, demog$rural), drop = TRUE)

expand_block <- function(b) {
  b <- b[match(old_ages, b$Age), ]            # canonical order
  pop_60_69 <- b$Population[b$Age == "60 to 69"]
  pop_70p   <- b$Population[b$Age == "70+"]
  # Splits (TODO: replace with ONS finer-grained values)
  pop_60_64 <- pop_60_69 * 0.5
  pop_65_69 <- pop_60_69 * 0.5
  pop_70_74 <- pop_70p   * (5 / 20)
  pop_75p   <- pop_70p   * (15 / 20)
  pop_65_74 <- pop_65_69 + pop_70_74


# Explicitly create 10 rows by duplicating row 9 as a scaffold for the new bands
new <- rbind(b, b[9, ])  # 10 rows; last 3 will be overwritten
  new$Age <- new_ages
  # 1..7 unchanged
  new$Population[1:7] <- b$Population[1:7]
  new$Population[8]   <- pop_60_64
  new$Population[9]   <- pop_65_74
  new$Population[10]  <- pop_75p
  new$IMD             <- b$IMD[1]
  new$rural           <- b$rural[1]
  new$tot_pop         <- sum(new$Population)
  new$Proportion      <- new$Population / new$tot_pop
  new
}

demog10 <- do.call(rbind, lapply(blocks, expand_block))
# Reorder: IMD asc, then rural (Urban first to match original), then age order
demog10 <- demog10[order(demog10$IMD,
                         factor(demog10$rural, levels = c("Urban","Rural")),
                         match(demog10$Age, new_ages)), ]
rownames(demog10) <- NULL

write.csv(demog10, file.path(output_dir, "demographics2021_10age.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s: %d rows (expected %d)\n",
            "data/demographics2021_10age.csv",
            nrow(demog10), 10 * length(unique(demog$IMD)) *
                          length(unique(demog$rural))))

## --- 2. Contact matrix -----------------------------------------------------
cm45 <- as.matrix(read.csv(file.path(input_dir, "Mas45_urban.csv"),
                           header = FALSE))
stopifnot(dim(cm45) == c(45, 45))

na_old <- 9
na_new <- 10
nimd   <- 5

# Indices 8=60-69, 9=70+; map old age idx -> new age idx (a vector of 1:10)
# new 1..7 = old 1..7 ; new 8 = old 8 ; new 9 = mean(old 8, old 9) ; new 10 = old 9
map_old_to_new <- function(mat9) {
  stopifnot(dim(mat9) == c(na_old, na_old))
  m <- matrix(0, na_new, na_new)
  # Row mapping: build expanded 10-row matrix from 9-row
  expand_rows <- rbind(mat9[1:7, , drop = FALSE],
                       mat9[8, , drop = FALSE],
                       (mat9[8, , drop = FALSE] + mat9[9, , drop = FALSE]) / 2,
                       mat9[9, , drop = FALSE])
  # Column mapping: same shape on the columns
  expand_cols <- cbind(expand_rows[, 1:7, drop = FALSE],
                       expand_rows[, 8, drop = FALSE],
                       (expand_rows[, 8, drop = FALSE] +
                        expand_rows[, 9, drop = FALSE]) / 2,
                       expand_rows[, 9, drop = FALSE])
  expand_cols
}

cm50 <- matrix(0, na_new * nimd, na_new * nimd)
for (ip in 1:nimd) {   # participant IMD block (rows)
  for (ic in 1:nimd) { # contact IMD block (cols)
    rows9 <- ((ip - 1) * na_old + 1):(ip * na_old)
    cols9 <- ((ic - 1) * na_old + 1):(ic * na_old)
    block <- cm45[rows9, cols9, drop = FALSE]
    rows10 <- ((ip - 1) * na_new + 1):(ip * na_new)
    cols10 <- ((ic - 1) * na_new + 1):(ic * na_new)
    cm50[rows10, cols10] <- map_old_to_new(block)
  }
}

write.table(cm50, file.path(output_dir, "Mas50_urban.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE)
cat(sprintf("Wrote %s: %d x %d\n",
            "data/Mas50_urban.csv", nrow(cm50), ncol(cm50)))
