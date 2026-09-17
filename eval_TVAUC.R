
source("R/setup.R")
source("R/functions_TVAUC.R")

# ----------------------------------------------------------------------------
# User settings
# ----------------------------------------------------------------------------

ntree_n <- 500
units <- 365.25
eval_years <- c(3, 6, 9, 12)
landmarkTimes <- eval_years * units

# ----------------------------------------------------------------------------
# Load test data and fitted models
# ----------------------------------------------------------------------------

load(file = paste0("Data/win_imp_sc_surv_", filename, ".Rdata"))

load(file = paste0("Data/cif_", filename, "_ntree", ntree_n, ".Rdata"))
cifmodel_fit <- cif_fit$ret

load(file = paste0("Data/lasso5cv_", filename, ".Rdata"))
load(file = paste0("Data/ridge5cv_", filename, ".Rdata"))
load(file = paste0("Data/elnet5cv_", filename, ".Rdata"))

load(file = paste0("Output/predictProb_result_", filename, ".Rdata"))

test_surv <- test_surv %>%
    arrange(ID, Start, Stop)

regularized_models <- list(
    LASSO = lasso_fit,
    Ridge = ridge_fit,
    Elastic_net = elnet_fit
)

# ----------------------------------------------------------------------------
# Obtain LTRC-CIF predictions
# ----------------------------------------------------------------------------
cif_risk_long <- extract_cif_risk_long(
    pred_object = pred_t,
    eval_times = landmarkTimes,
    eval_labels = eval_years
)

# ----------------------------------------------------------------------------
# Build one combined score data frame at each evaluation year
# ----------------------------------------------------------------------------
combined_score_list <- purrr::map(eval_years, function(yr) {
    t0 <- yr * units

    # Active covariate row for each subject at t0
    lm_rows <- make_landmark_rows(test_surv, t0)
    active_ids <- as.character(lm_rows$ID)

    # Regularized Cox linear predictors
    reg_scores <- tibble(
        ID = active_ids,
        year = yr,
        eval_time = t0
    )

    for (model_name in names(regularized_models)) {
        reg_scores[[model_name]] <- as.numeric(
            predict(
                regularized_models[[model_name]],
                newdata = lm_rows,
                type = "lp"
            )
        )
    }

    # LTRC-CIF scores for subjects active at t0
    cif_scores <- cif_risk_long %>%
        filter(
            year == yr,
            as.character(ID) %in% active_ids
        ) %>%
        dplyr::select(ID, LTRC_CIF)

    inner_join(reg_scores, cif_scores, by = "ID")
})

names(combined_score_list) <- paste0(eval_years, "_year")

# ----------------------------------------------------------------------------
# Calculate time-dependent AUC
# ----------------------------------------------------------------------------

tdcauc_check <- get_tdcauc(
    combined_score_list = combined_score_list,
    survival_data = test_surv,
    eval_years = eval_years,
    units = units,
    method = "Cox"
)

print(tdcauc_check$AUC_wide)

save(
    tdcauc_check,
    file = "Output/tdcauc.Rdata"
)

# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
p <- ggplot(
    tdcauc_check$AUC,
    aes(
        x = times,
        y = estimate,
        colour = model,
        group = model
    )
) +
    geom_line() +
    geom_point() +
    labs(
        x = "Evaluation time (years)",
        y = "Time-dependent AUC",
        colour = "Model"
    ) +
    theme_minimal()

p
