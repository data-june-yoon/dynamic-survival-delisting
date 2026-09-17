
# ============================================================================
# Time-dependent incident/dynamic AUC at 3, 6, 9, and 12 years
# Models: LTRC-CIF, LASSO Cox, Ridge Cox, Elastic-net Cox
#
# Main design:
#   1. Use one active covariate record per subject at each evaluation time.
#   2. Use Cox linear predictors from the interval active at time t.
#   3. Transform LTRC-CIF as log(-log(S_hat(t))) for Cox-based risksetROC.
#   4. Compare all models on the same subjects at each evaluation time.
#   5. Retain each subject's original delayed-entry time for risksetROC.
# ============================================================================


extract_cif_risk_long <- function(pred_object, eval_times, eval_labels) {
    
    surv_time_by_subject <- pred_object$survival.probs

    # survival.id may contain the same subject ID repeatedly across
    # prediction times. Retain one ID per subject, as in v2 and v3.
    ids <- unique(as.character(pred_object$survival.id))

    if (any(surv_time_by_subject < 0 | surv_time_by_subject > 1,
            na.rm = TRUE)) {
        warning("Survival probabilities out of range [0, 1]")
    }

    risk_time_by_subject <- 1 - surv_time_by_subject

    colnames(risk_time_by_subject) <- ids
    rownames(risk_time_by_subject) <- as.character(eval_labels)

    as.data.frame(t(risk_time_by_subject), check.names = FALSE) %>%
        rownames_to_column("ID") %>%
        pivot_longer(
            cols = -ID,
            names_to = "year",
            values_to = "LTRC_CIF"
        ) %>%
        mutate(
            ID = as.character(ID),
            year = as.numeric(year)
        )
}

# ----------------------------------------------------------------------------
# Return one covariate row per subject that is active at t0.
# Counting-process convention: Start < t0 <= Stop.
# ----------------------------------------------------------------------------
make_landmark_rows <- function(data, t0) {
    data %>%
        filter(Start < t0, Stop >= t0) %>%
        group_by(ID) %>%
        arrange(Start, Stop, .by_group = TRUE) %>%
        slice_tail(n = 1L) %>%
        ungroup()
}

# ----------------------------------------------------------------------------
# Create one subject-level entry, follow-up time, and terminal event status.
# The subject's original entry time is retained rather than the active
# time-varying interval start.
# ----------------------------------------------------------------------------
make_subject_outcomes <- function(data,
                                  idvar = "ID",
                                  entryvar = "Start",
                                  timevar = "Stop",
                                  statusvar = "Event") {
    data %>%
        group_by(.data[[idvar]]) %>%
        summarise(
            entry = min(.data[[entryvar]], na.rm = TRUE),
            ultimate_time = max(.data[[timevar]], na.rm = TRUE),
            status = as.integer(any(.data[[statusvar]] == 1L,
                                    na.rm = TRUE)),
            .groups = "drop"
        ) %>%
        mutate(ID = as.character(.data[[idvar]]))
}

# ----------------------------------------------------------------------------
# AUC for one pre-computed marker at one evaluation time.
# ----------------------------------------------------------------------------
one_tdcauc <- function(data, marker_col, predict.time,
                       method = "Cox", ...) {
    dat <- data %>%
        filter(
            is.finite(entry),
            is.finite(ultimate_time),
            is.finite(.data[[marker_col]]),
            status %in% c(0L, 1L),
            entry < predict.time,
            ultimate_time >= predict.time
        )

    if (nrow(dat) < 2L ||
        length(unique(dat[[marker_col]])) < 2L ||
        sum(dat$status == 1L, na.rm = TRUE) == 0L) {
        return(NA_real_)
    }

    auc_fit <- risksetROC::risksetROC(
        Stime = dat$ultimate_time,
        status = dat$status,
        entry = dat$entry,
        marker = dat[[marker_col]],
        predict.time = predict.time,
        method = method,
        plot = FALSE,
        ...
    )

    as.numeric(auc_fit$AUC[1L])
}

# ----------------------------------------------------------------------------
# AUCs for all evaluation years and marker columns.
# ----------------------------------------------------------------------------
get_tdcauc <- function(combined_score_list,
                       survival_data,
                       eval_years,
                       units = 365.25,
                       method = "Cox",
                       marker_cols = c(
                           "LTRC_CIF",
                           "LASSO",
                           "Ridge",
                           "Elastic_net"
                       ),
                       ...) {

    subject_outcomes <- make_subject_outcomes(survival_data)

    all_scores <- bind_rows(combined_score_list) %>%
        distinct(ID, year, .keep_all = TRUE)

    auc_data <- inner_join(subject_outcomes, all_scores, by = "ID")

    results <- list()

    for (yr in eval_years) {
        t0 <- yr * units

        risk_set <- auc_data %>%
            filter(
                year == yr,
                entry < t0,
                ultimate_time >= t0
            )

        # Transform the LTRC-CIF marker for Cox-based risksetROC.
        if (method == "Cox") {
            eps <- 1e-8
            risk_set <- risk_set %>%
                mutate(
                    LTRC_CIF = pmin(pmax(LTRC_CIF, eps), 1 - eps),
                    LTRC_CIF = log(-log(1 - LTRC_CIF))
                )
        }

        # All models are evaluated using the same subjects at this year.
        risk_set <- risk_set %>%
            filter(if_all(all_of(marker_cols), is.finite))

        for (model_name in marker_cols) {
            auc_est <- one_tdcauc(
                data = risk_set,
                marker_col = model_name,
                predict.time = t0,
                method = method,
                ...
            )

            results[[length(results) + 1L]] <- tibble(
                model = model_name,
                times = yr,
                year = yr,
                eval_time = t0,
                estimate = auc_est,
                n = nrow(risk_set),
                events = sum(risk_set$status == 1L, na.rm = TRUE)
            )
        }
    }

    AUC <- bind_rows(results) %>%
        mutate(model = factor(model, levels = marker_cols)) %>%
        arrange(year, model) %>%
        mutate(model = as.character(model))

    AUC_wide <- AUC %>%
        dplyr::select(model, year, estimate) %>%
        pivot_wider(
            names_from = year,
            values_from = estimate,
            names_prefix = "year_"
        ) %>%
        arrange(match(model, marker_cols))

    out <- list(
        AUC = as.data.frame(AUC),
        AUC_wide = as.data.frame(AUC_wide)
    )

    class(out) <- "tdcauc"
    return(out)
}
