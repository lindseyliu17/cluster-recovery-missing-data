# =============================================================================
# main.R
#
# Entry point for the simulation study:
#   "Which Proximity Measures Are Most Robust to Missing Data
#    for Recovering Cluster Structure?"
#
# HOW TO RUN:
#   Open this file in RStudio and click "Source", or run in the terminal:
#     Rscript main.R
#   Make sure your working directory is the project root (the folder
#   containing this file).
#
# WHAT THIS SCRIPT DOES:
#   0. Installs / loads required packages
#   1. Sources all module files from the R/ folder
#   2. Defines simulation conditions
#   3. Runs the full simulation (saves results to output/)
#   4. Generates all tables and figures
#
# RUNTIME NOTE:
#   With N_SIM = 300 replicates × 28 conditions = 8,400 replications,
#   each running 10 clustering operations, expect roughly 30–90 minutes
#   or longer on a typical laptop. To do a quick test first, set
#   RUN_DRY_RUN <- TRUE in Section 3a.
# =============================================================================


# =============================================================================
# 0. Package setup
# =============================================================================

required_packages <- c(
  "MASS",       # mvrnorm() — multivariate normal simulation
  "cluster",    # pam(), daisy(), silhouette()
  "mclust",     # adjustedRandIndex()
  "ggplot2",    # all plotting
  "dplyr",      # data frame manipulation
  "tidyr",      # pivoting / reshaping
  "scales"      # percent labels on ggplot axes
)

# Install any missing packages automatically
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

cat("All packages loaded.\n\n")


# =============================================================================
# 1. Source all module files
# =============================================================================

source("R/01_generate_data.R")   # generate_cluster_data()
source("R/02_missingness.R")     # introduce_missingness()
source("R/03_proximity.R")       # compute_all_proximity()
source("R/04_clustering.R")      # cluster_from_distance()
source("R/05_evaluation.R")      # compute_ari(), compute_silhouette()
source("R/06_simulation.R")      # run_one_replicate(), run_full_simulation()
source("R/07_results.R")         # create_all_outputs()

cat("All modules sourced.\n\n")


# =============================================================================
# 2. Create output directories (safe to re-run: showWarnings = FALSE)
# =============================================================================

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 3. Define simulation conditions
# =============================================================================

# ---- Fixed settings ----------------------------------------------------------
N_PER_CLUSTER <- 100    # observations per cluster → total n = 300
N_VARS        <- 5      # number of clustering variables in each dataset
N_CLUSTERS    <- 3      # number of true clusters (also passed to algorithms)
N_SIM         <- 300    # replicates per condition (reduce for testing)
SEED          <- 5124   # global random seed for reproducibility

# ---- Complete-data baseline (no missingness) ---------------------------------
# Including a baseline lets us see full-information performance for reference.
baseline_conditions <- data.frame(
  data_type  = rep(c("continuous", "mixed"), each = 2),
  separation = rep(c("moderate", "strong"), times = 2),
  mechanism  = "none",
  miss_rate  = 0,
  stringsAsFactors = FALSE
)

# ---- Missing-data conditions -------------------------------------------------
# Full factorial:
#   2 data types × 2 separations × 2 mechanisms × 3 miss rates = 24 conditions
missing_conditions <- expand.grid(
  data_type  = c("continuous", "mixed"),
  separation = c("moderate", "strong"),
  mechanism  = c("MCAR", "MAR"),
  miss_rate  = c(0.10, 0.30, 0.50),
  stringsAsFactors = FALSE
)

# Combine into one conditions table (28 rows total)
sim_conditions <- rbind(baseline_conditions, missing_conditions)

cat("Simulation design summary:\n")
cat(sprintf("  Conditions       : %d\n", nrow(sim_conditions)))
cat(sprintf("  Replicates/cond  : %d\n", N_SIM))
cat(sprintf("  Total replications: %d\n", nrow(sim_conditions) * N_SIM))
cat(sprintf("  Data types        : 2  (continuous, mixed)\n"))
cat(sprintf("  Proximity measures: 5\n"))
cat(sprintf("  Clustering methods: 2  (PAM, UPGMA / Average Linkage)\n"))
cat(sprintf("  Evaluations total : %d\n\n", nrow(sim_conditions) * N_SIM * 5 * 2))


# =============================================================================
# 3a. (OPTIONAL) Dry run — quick check that everything runs without errors
#     Set RUN_DRY_RUN <- TRUE to execute; FALSE skips it and goes straight
#     to the full simulation.
# =============================================================================

RUN_DRY_RUN <- FALSE   # ← change to TRUE to run the quick test first

if (RUN_DRY_RUN) {
  cat("--- DRY RUN (5 replicates per condition) ---\n")
  set.seed(SEED)
  test_results <- run_full_simulation(
    conditions    = sim_conditions,
    n_sim         = 5,           # tiny number just to verify flow
    n_per_cluster = N_PER_CLUSTER,
    n_vars        = N_VARS,
    n_clusters    = N_CLUSTERS
  )
  cat("Dry run complete. Rows in results:", nrow(test_results), "\n\n")
}


# =============================================================================
# 4. Run full simulation
# =============================================================================

cat("Starting full simulation...\n")
start_time <- proc.time()

set.seed(SEED)
results <- run_full_simulation(
  conditions    = sim_conditions,
  n_sim         = N_SIM,
  n_per_cluster = N_PER_CLUSTER,
  n_vars        = N_VARS,
  n_clusters    = N_CLUSTERS
)

elapsed <- proc.time() - start_time
cat(sprintf("\nTotal runtime: %.1f minutes\n\n", elapsed["elapsed"] / 60))

# Save raw results so you can reload them later without re-running
saveRDS(results, "output/simulation_results.rds")
cat("Raw results saved to: output/simulation_results.rds\n\n")


# =============================================================================
# 5. Generate tables and figures
# =============================================================================
# To re-run just this section after reloading saved results:
#   results <- readRDS("output/simulation_results.rds")
#   create_all_outputs(results)

cat("Generating tables and figures...\n")
create_all_outputs(results)

cat("\nDone! Check output/figures/ and output/tables/\n")
