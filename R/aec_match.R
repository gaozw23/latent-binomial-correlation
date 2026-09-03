haversine_km <- function(lat1, lon1, lat2, lon2) {
  values <- c(lat1, lon1, lat2, lon2)
  if (length(values) != 4L || any(!is.finite(values))) stop("Coordinates must be four finite scalars", call. = FALSE)
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

match_aec_places <- function(...) {
  audit_aec_matching_feasibility(...)
}

.aec_place_population <- function(metadata) {
  x <- metadata[metadata$polling_place_type_id == "1", , drop = FALSE]
  x$latitude[!x$coordinate_valid] <- NA_real_
  x$longitude[!x$coordinate_valid] <- NA_real_
  x$name_normalised <- normalise_place_name(x$polling_place_name)
  x$key_id <- paste(x$state_ab, x$polling_place_id, sep = "|")
  x$key_name <- paste(x$state_ab, x$name_normalised, sep = "|")
  id_frequency <- table(x$key_id)
  x$id_unique <- unname(id_frequency[x$key_id]) == 1L
  x
}

.aec_distance <- function(lat1, lon1, lat2, lon2) {
  valid <- is.finite(lat1) & is.finite(lon1) & is.finite(lat2) & is.finite(lon2)
  out <- rep(NA_real_, length(lat1))
  if (any(valid)) {
    rad <- pi / 180
    dlat <- (lat2[valid] - lat1[valid]) * rad
    dlon <- (lon2[valid] - lon1[valid]) * rad
    a <- sin(dlat / 2)^2 + cos(lat1[valid] * rad) * cos(lat2[valid] * rad) * sin(dlon / 2)^2
    out[valid] <- 6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
  }
  out
}

audit_aec_matching_feasibility <- function(audit) {
  p19 <- .aec_place_population(audit$inventory$metadata[["2019"]])
  p22 <- .aec_place_population(audit$inventory$metadata[["2022"]])
  states <- c(.aec_jurisdictions, "NATIONAL")
  feasibility <- list()
  exact_candidates <- list()
  secondary_candidates <- list()
  ambiguous_rows <- list()
  fz <- ez <- sz <- az <- 0L

  for (state in states) {
    a <- if (state == "NATIONAL") p19 else p19[p19$state_ab == state, , drop = FALSE]
    b <- if (state == "NATIONAL") p22 else p22[p22$state_ab == state, , drop = FALSE]
    a_unique <- a[a$id_unique, , drop = FALSE]
    b_unique <- b[b$id_unique, , drop = FALSE]
    exact_keys <- intersect(a_unique$key_id, b_unique$key_id)
    exact <- merge(
      a_unique[a_unique$key_id %in% exact_keys, , drop = FALSE],
      b_unique[b_unique$key_id %in% exact_keys, , drop = FALSE],
      by = "key_id", suffixes = c("_2019", "_2022"), sort = FALSE
    )
    if (nrow(exact)) {
      exact$distance_km <- .aec_distance(exact$latitude_2019, exact$longitude_2019,
                                         exact$latitude_2022, exact$longitude_2022)
      exact$coordinate_support <- !is.na(exact$distance_km)
      exact$coordinate_over_5km <- exact$coordinate_support & exact$distance_km > 5
    } else {
      exact$distance_km <- numeric()
      exact$coordinate_support <- logical()
      exact$coordinate_over_5km <- logical()
    }
    if (state != "NATIONAL") {
      ez <- ez + 1L
      exact_candidates[[ez]] <- exact[c(
        "state_ab_2019", "polling_place_id_2019", "polling_place_name_2019",
        "polling_place_name_2022", "distance_km", "coordinate_support", "coordinate_over_5km"
      )]
    }

    unmatched_a <- a[!a$key_id %in% exact_keys, , drop = FALSE]
    unmatched_b <- b[!b$key_id %in% exact_keys, , drop = FALSE]
    name_count_a <- table(unmatched_a$key_name)
    name_count_b <- table(unmatched_b$key_name)
    overlap_names <- intersect(names(name_count_a), names(name_count_b))
    unique_names <- overlap_names[name_count_a[overlap_names] == 1L & name_count_b[overlap_names] == 1L]
    ambiguous_names <- setdiff(overlap_names, unique_names)
    secondary <- merge(
      unmatched_a[unmatched_a$key_name %in% unique_names, , drop = FALSE],
      unmatched_b[unmatched_b$key_name %in% unique_names, , drop = FALSE],
      by = "key_name", suffixes = c("_2019", "_2022"), sort = FALSE
    )
    if (nrow(secondary)) {
      secondary$distance_km <- .aec_distance(
        secondary$latitude_2019, secondary$longitude_2019,
        secondary$latitude_2022, secondary$longitude_2022
      )
      secondary$coordinate_support <- !is.na(secondary$distance_km)
      secondary$coordinate_pass_2km <- secondary$coordinate_support & secondary$distance_km <= 2
      secondary$coordinate_fail_2km <- secondary$coordinate_support & secondary$distance_km > 2
      secondary$potentially_acceptable <- secondary$coordinate_pass_2km | !secondary$coordinate_support
    } else {
      secondary$distance_km <- numeric()
      secondary$coordinate_support <- logical()
      secondary$coordinate_pass_2km <- logical()
      secondary$coordinate_fail_2km <- logical()
      secondary$potentially_acceptable <- logical()
    }
    if (state != "NATIONAL") {
      sz <- sz + 1L
      secondary_candidates[[sz]] <- secondary[c(
        "state_ab_2019", "polling_place_id_2019", "polling_place_name_2019",
        "polling_place_id_2022", "polling_place_name_2022", "name_normalised_2019",
        "latitude_2019", "longitude_2019", "latitude_2022", "longitude_2022",
        "distance_km", "coordinate_support", "coordinate_pass_2km",
        "coordinate_fail_2km", "potentially_acceptable"
      )]
      if (length(ambiguous_names)) {
        az <- az + 1L
        ambiguous_rows[[az]] <- data.frame(
          state_ab = state,
          name_normalised = sub("^[^|]+[|]", "", ambiguous_names),
          count_2019 = as.integer(name_count_a[ambiguous_names]),
          count_2022 = as.integer(name_count_b[ambiguous_names]),
          stringsAsFactors = FALSE
        )
      }
    }

    duplicated_candidates <- sum(duplicated(secondary$polling_place_id_2019)) +
      sum(duplicated(secondary$polling_place_id_2022))
    fz <- fz + 1L
    feasibility[[fz]] <- data.frame(
      state_ab = state, ordinary_places_2019 = nrow(a), ordinary_places_2022 = nrow(b),
      exact_id_match_candidates = nrow(exact),
      exact_id_percent_of_2019 = if (nrow(a)) 100 * nrow(exact) / nrow(a) else NA_real_,
      exact_id_over_5km = sum(exact$coordinate_over_5km, na.rm = TRUE),
      exact_id_coordinate_missing = sum(!exact$coordinate_support),
      unmatched_2019_after_exact_id = nrow(unmatched_a),
      unmatched_2022_after_exact_id = nrow(unmatched_b),
      potential_unique_name_matches = nrow(secondary),
      secondary_coordinate_evaluated = sum(secondary$coordinate_support),
      secondary_coordinate_pass_2km = sum(secondary$coordinate_pass_2km),
      secondary_coordinate_fail_2km = sum(secondary$coordinate_fail_2km),
      secondary_coordinate_missing = sum(!secondary$coordinate_support),
      potential_secondary_acceptable = sum(secondary$potentially_acceptable),
      ambiguous_name_keys = length(ambiguous_names),
      duplicated_candidate_matches = duplicated_candidates,
      remaining_unmatched_2019_after_potential_secondary = nrow(unmatched_a) - sum(secondary$potentially_acceptable),
      remaining_unmatched_2022_after_potential_secondary = nrow(unmatched_b) - sum(secondary$potentially_acceptable),
      stringsAsFactors = FALSE
    )
  }

  secondary <- if (length(secondary_candidates)) do.call(rbind, secondary_candidates) else data.frame()
  exact <- if (length(exact_candidates)) do.call(rbind, exact_candidates) else data.frame()
  ambiguous <- if (length(ambiguous_rows)) do.call(rbind, ambiguous_rows) else
    data.frame(state_ab = character(), name_normalised = character(),
               count_2019 = integer(), count_2022 = integer())
  distance <- secondary$distance_km[is.finite(secondary$distance_km)]
  probabilities <- c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
  distance_distribution <- data.frame(
    probability = probabilities,
    distance_km = if (length(distance)) as.numeric(quantile(distance, probabilities, type = 7)) else NA_real_,
    stringsAsFactors = FALSE
  )
  near_threshold <- secondary[is.finite(secondary$distance_km) &
                                secondary$distance_km >= 1.5 & secondary$distance_km <= 2.5, , drop = FALSE]
  state_feasibility <- do.call(rbind, feasibility)
  matching_anomalies <- do.call(rbind, lapply(seq_len(nrow(state_feasibility) - 1L), function(i) {
    state <- state_feasibility$state_ab[i]
    rbind(
      .aec_anomaly_row(NA_integer_, state, "matching_feasibility", "exact_id_coordinate_over_5km",
                       state_feasibility$exact_id_over_5km[i],
                       "Exact-ID candidate retained only for manual audit; not accepted in Phase 5A.",
                       if (state_feasibility$exact_id_over_5km[i] > 0) "REVIEW" else "PASS"),
      .aec_anomaly_row(NA_integer_, state, "matching_feasibility", "secondary_coordinate_over_2km",
                       state_feasibility$secondary_coordinate_fail_2km[i],
                       "Unique-name candidate fails the 2 km rule and is not potentially acceptable.",
                       if (state_feasibility$secondary_coordinate_fail_2km[i] > 0) "REVIEW" else "PASS"),
      .aec_anomaly_row(NA_integer_, state, "matching_feasibility", "ambiguous_normalised_name_key",
                       state_feasibility$ambiguous_name_keys[i],
                       "Non-unique normalised-name keys are not candidates for automatic matching.",
                       if (state_feasibility$ambiguous_name_keys[i] > 0) "REVIEW" else "PASS")
    )
  }))
  list(feasibility = do.call(rbind, feasibility), exact_candidates = exact,
       secondary_candidates = secondary, ambiguous_names = ambiguous,
       distance_distribution = distance_distribution, near_threshold = near_threshold,
       anomaly_inventory = matching_anomalies)
}

.aec_full_address <- function(x) {
  fields <- c("premises_name", "address_1", "address_2", "address_3",
              "suburb", "premises_state_ab", "postcode")
  apply(x[fields], 1L, function(parts) paste(parts[nzchar(parts)], collapse = ", "))
}

.aec_token_jaccard <- function(a, b) {
  tokens <- function(x) unique(strsplit(normalise_place_name(x), " ", fixed = TRUE)[[1L]])
  aa <- tokens(a)
  bb <- tokens(b)
  if (!length(union(aa, bb))) return(1)
  length(intersect(aa, bb)) / length(union(aa, bb))
}

build_aec_exact_id_outlier_adjudication <- function(audit, disposition_path = NULL) {
  p19 <- .aec_place_population(audit$inventory$metadata[["2019"]])
  p22 <- .aec_place_population(audit$inventory$metadata[["2022"]])
  exact <- merge(p19, p22, by = "key_id", suffixes = c("_2019", "_2022"), sort = FALSE)
  exact$distance_km <- .aec_distance(exact$latitude_2019, exact$longitude_2019,
                                      exact$latitude_2022, exact$longitude_2022)
  zero_coordinate_ids <- c("83845", "2887")
  selected <- exact[(!is.na(exact$distance_km) & exact$distance_km > 5) |
                      exact$polling_place_id_2022 %in% zero_coordinate_ids, , drop = FALSE]
  if (nrow(selected) != 27L) stop("Expected exactly 27 exact-ID adjudication records", call. = FALSE)
  address_fields <- c("premises_name", "address_1", "address_2", "address_3",
                      "suburb", "premises_state_ab", "postcode")
  assemble <- function(i, suffix) {
    values <- vapply(address_fields, function(field) selected[[paste0(field, suffix)]][i], character(1))
    paste(values[nzchar(values)], collapse = ", ")
  }
  selected$full_address_2019 <- vapply(seq_len(nrow(selected)), assemble, character(1), suffix = "_2019")
  selected$full_address_2022 <- vapply(seq_len(nrow(selected)), assemble, character(1), suffix = "_2022")
  selected$normalised_name_equal <- selected$name_normalised_2019 == selected$name_normalised_2022
  selected$exact_original_name_equal <- selected$polling_place_name_2019 == selected$polling_place_name_2022
  selected$normalised_premises_equal <- normalise_place_name(selected$premises_name_2019) ==
    normalise_place_name(selected$premises_name_2022)
  selected$normalised_address_equal <- normalise_place_name(selected$full_address_2019) ==
    normalise_place_name(selected$full_address_2022)
  selected$address_token_jaccard <- vapply(seq_len(nrow(selected)), function(i) {
    .aec_token_jaccard(selected$full_address_2019[i], selected$full_address_2022[i])
  }, numeric(1))
  selected$relevant_address_evidence <- sprintf(
    "same premises=%s; same normalized full address=%s; address-token Jaccard=%.3f; same division=%s",
    selected$normalised_premises_equal, selected$normalised_address_equal,
    selected$address_token_jaccard, selected$division_id_2019 == selected$division_id_2022
  )
  selected$id_unique_2019 <- selected$id_unique_2019
  selected$id_unique_2022 <- selected$id_unique_2022
  selected$proposed_disposition <- "UNRESOLVED"
  selected$reason <- "Awaiting structured adjudication"
  if (!is.null(disposition_path)) {
    disposition <- read.csv(disposition_path, stringsAsFactors = FALSE, colClasses = "character")
    required <- c("state_ab", "polling_place_id", "proposed_disposition", "reason")
    if (!all(required %in% names(disposition)) || nrow(disposition) != 27L ||
        anyDuplicated(disposition[c("state_ab", "polling_place_id")])) {
      stop("Outlier disposition file must contain exactly 27 unique records", call. = FALSE)
    }
    key <- paste(selected$state_ab_2019, selected$polling_place_id_2019, sep = "|")
    disposition_key <- paste(disposition$state_ab, disposition$polling_place_id, sep = "|")
    order <- match(key, disposition_key)
    if (anyNA(order)) stop("A required outlier disposition is missing", call. = FALSE)
    selected$proposed_disposition <- disposition$proposed_disposition[order]
    selected$reason <- disposition$reason[order]
  }
  allowed <- c("ACCEPT_EXACT_ID", "REJECT_IDENTITY_CHANGE", "UNRESOLVED")
  if (any(!selected$proposed_disposition %in% allowed)) stop("Invalid outlier disposition", call. = FALSE)
  out <- data.frame(
    State = selected$state_ab_2019,
    PollingPlaceID = selected$polling_place_id_2019,
    DivisionID_2019 = selected$division_id_2019,
    DivisionNm_2019 = selected$division_name_2019,
    DivisionID_2022 = selected$division_id_2022,
    DivisionNm_2022 = selected$division_name_2022,
    PollingPlaceNm_2019 = selected$polling_place_name_2019,
    PollingPlaceNm_2022 = selected$polling_place_name_2022,
    PremisesNm_2019 = selected$premises_name_2019,
    PremisesNm_2022 = selected$premises_name_2022,
    FullAddress_2019 = selected$full_address_2019,
    FullAddress_2022 = selected$full_address_2022,
    Latitude_2019 = selected$latitude_text_2019,
    Longitude_2019 = selected$longitude_text_2019,
    Latitude_2022 = selected$latitude_text_2022,
    Longitude_2022 = selected$longitude_text_2022,
    GreatCircleDistanceKm = selected$distance_km,
    NormalizedNameEquality = selected$normalised_name_equal,
    ExactOriginalNameEquality = selected$exact_original_name_equal,
    RelevantAddressEvidence = selected$relevant_address_evidence,
    PollingPlaceIDUnique2019 = selected$id_unique_2019,
    PollingPlaceIDUnique2022 = selected$id_unique_2022,
    ProposedDisposition = selected$proposed_disposition,
    Reason = selected$reason,
    stringsAsFactors = FALSE
  )
  out <- out[order(match(out$State, .aec_jurisdictions), as.numeric(out$PollingPlaceID)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

aggregate_aec_polling_place_year <- function(audit, year) {
  year <- as.integer(year)
  result_keys <- startsWith(names(audit$inventory$results), paste0(year, "_"))
  results <- do.call(rbind, audit$inventory$results[result_keys])
  metadata <- audit$inventory$metadata[[as.character(year)]]
  metadata <- metadata[metadata$polling_place_type_id == "1", , drop = FALSE]
  metadata$key_id <- paste(metadata$state_ab, metadata$polling_place_id, sep = "|")
  if (anyDuplicated(metadata$key_id)) stop("Ordinary metadata state/ID is not unique", call. = FALSE)
  results$key_id <- paste(results$state_ab, results$polling_place_id, sep = "|")
  results <- results[results$key_id %in% metadata$key_id, , drop = FALSE]
  split_rows <- split(seq_len(nrow(results)), results$key_id)
  rows <- lapply(split_rows, function(ii) {
    x <- results[ii, , drop = FALSE]
    informal <- !nzchar(x$party_ab) & toupper(x$party_nm) == "INFORMAL"
    formal <- !informal
    if (sum(informal) != 1L) stop("Expected exactly one informal row per ordinary place", call. = FALSE)
    alp <- formal & toupper(x$party_ab) == "ALP"
    if (sum(alp) != 1L) stop("Expected exactly one ALP row per ordinary place", call. = FALSE)
    data.frame(
      key_id = x$key_id[1L], state_ab = x$state_ab[1L],
      division_id_results = x$division_id[1L], division_name_results = x$division_name[1L],
      polling_place_id = x$polling_place_id[1L], polling_place_name_results = x$polling_place_name[1L],
      X = sum(x$ordinary_votes[alp]), M = sum(x$ordinary_votes[formal]),
      stringsAsFactors = FALSE
    )
  })
  aggregated <- do.call(rbind, rows)
  rownames(aggregated) <- NULL
  keep <- c("key_id", "state_ab", "division_id", "division_name", "polling_place_id",
            "polling_place_name", "name_normalised", "latitude", "longitude",
            "coordinate_valid", "premises_name", "address_1", "address_2", "address_3",
            "suburb", "premises_state_ab", "postcode")
  metadata$name_normalised <- normalise_place_name(metadata$polling_place_name)
  out <- merge(aggregated, metadata[keep], by = c("key_id", "state_ab", "polling_place_id"),
               all.x = TRUE, sort = FALSE)
  if (nrow(out) != nrow(metadata) || anyNA(out$division_id)) {
    stop("Ordinary result-to-metadata aggregation is incomplete", call. = FALSE)
  }
  if (any(!is.finite(out$X)) || any(!is.finite(out$M)) ||
      any(out$X < 0 | out$X > out$M) || any(out$M != round(out$M))) {
    stop("Aggregated AEC count constraints failed", call. = FALSE)
  }
  out$Y <- ifelse(out$M > 0, out$X / out$M, NA_real_)
  out
}

freeze_aec_matched_table <- function(audit, adjudication,
                                     frozen_path = "results/derived/aec_matched_frozen.csv",
                                     exclusion_path = "results/derived/aec_match_exclusions.csv") {
  if (any(adjudication$ProposedDisposition == "UNRESOLVED")) {
    stop("Cannot freeze matches while an exact-ID adjudication remains unresolved", call. = FALSE)
  }
  year19 <- aggregate_aec_polling_place_year(audit, 2019L)
  year22 <- aggregate_aec_polling_place_year(audit, 2022L)
  common <- intersect(year19$key_id, year22$key_id)
  rejected <- adjudication[adjudication$ProposedDisposition == "REJECT_IDENTITY_CHANGE", , drop = FALSE]
  rejected_keys <- paste(rejected$State, rejected$PollingPlaceID, sep = "|")
  accepted_exact_keys <- setdiff(common, rejected_keys)

  exact <- data.frame(
    match_route = "EXACT_ID",
    state_ab = sub("[|].*$", "", accepted_exact_keys),
    polling_place_id_2019 = sub("^[^|]+[|]", "", accepted_exact_keys),
    polling_place_id_2022 = sub("^[^|]+[|]", "", accepted_exact_keys),
    stringsAsFactors = FALSE
  )
  matching <- audit_aec_matching_feasibility(audit)
  secondary_source <- matching$secondary_candidates
  secondary_source <- secondary_source[secondary_source$potentially_acceptable, , drop = FALSE]
  if (nrow(secondary_source) != 15L || any(!secondary_source$coordinate_pass_2km)) {
    stop("Expected exactly 15 coordinate-supported secondary candidates", call. = FALSE)
  }
  secondary <- data.frame(
    match_route = "UNIQUE_NAME_COORD", state_ab = secondary_source$state_ab_2019,
    polling_place_id_2019 = secondary_source$polling_place_id_2019,
    polling_place_id_2022 = secondary_source$polling_place_id_2022,
    stringsAsFactors = FALSE
  )
  candidates <- rbind(exact, secondary)
  if (anyDuplicated(paste(candidates$state_ab, candidates$polling_place_id_2019)) ||
      anyDuplicated(paste(candidates$state_ab, candidates$polling_place_id_2022))) {
    stop("Candidate matched table is not one-to-one", call. = FALSE)
  }

  names(year19) <- paste0(names(year19), "_2019")
  names(year22) <- paste0(names(year22), "_2022")
  joined <- merge(candidates, year19,
                  by.x = c("state_ab", "polling_place_id_2019"),
                  by.y = c("state_ab_2019", "polling_place_id_2019"), all.x = TRUE, sort = FALSE)
  joined <- merge(joined, year22,
                  by.x = c("state_ab", "polling_place_id_2022"),
                  by.y = c("state_ab_2022", "polling_place_id_2022"), all.x = TRUE, sort = FALSE)
  if (nrow(joined) != nrow(candidates) || anyNA(joined$M_2019) || anyNA(joined$M_2022)) {
    stop("Failed to attach year-specific aggregates to match candidates", call. = FALSE)
  }
  joined$match_distance_km <- .aec_distance(
    ifelse(joined$coordinate_valid_2019, joined$latitude_2019, NA_real_),
    ifelse(joined$coordinate_valid_2019, joined$longitude_2019, NA_real_),
    ifelse(joined$coordinate_valid_2022, joined$latitude_2022, NA_real_),
    ifelse(joined$coordinate_valid_2022, joined$longitude_2022, NA_real_)
  )
  adjudication_key <- paste(adjudication$State, adjudication$PollingPlaceID, sep = "|")
  joined_key <- paste(joined$state_ab, joined$polling_place_id_2019, sep = "|")
  adjudication_index <- match(joined_key, adjudication_key)
  joined$exact_id_outlier_adjudication <- ifelse(
    joined$match_route == "UNIQUE_NAME_COORD", "NOT_APPLICABLE",
    ifelse(is.na(adjudication_index), "NOT_REQUIRED_DISTANCE_WITHIN_5KM",
           adjudication$ProposedDisposition[adjudication_index])
  )
  joined$eligible <- joined$M_2019 >= 2 & joined$M_2022 >= 2

  identity_exclusions <- data.frame(
    State = rejected$State, PollingPlaceID_2019 = rejected$PollingPlaceID,
    PollingPlaceID_2022 = rejected$PollingPlaceID, MatchRoute = "EXACT_ID",
    M_2019 = NA_real_, M_2022 = NA_real_,
    ExclusionReason = paste0("REJECT_IDENTITY_CHANGE: ", rejected$Reason),
    stringsAsFactors = FALSE
  )
  denominator_exclusions <- data.frame(
    State = joined$state_ab[!joined$eligible],
    PollingPlaceID_2019 = joined$polling_place_id_2019[!joined$eligible],
    PollingPlaceID_2022 = joined$polling_place_id_2022[!joined$eligible],
    MatchRoute = joined$match_route[!joined$eligible],
    M_2019 = joined$M_2019[!joined$eligible], M_2022 = joined$M_2022[!joined$eligible],
    ExclusionReason = "M_LT_2_IN_AT_LEAST_ONE_YEAR",
    stringsAsFactors = FALSE
  )
  exclusions <- rbind(identity_exclusions, denominator_exclusions)

  frozen <- joined[joined$eligible, , drop = FALSE]
  frozen <- frozen[order(match(frozen$state_ab, .aec_jurisdictions),
                         frozen$polling_place_id_2019, frozen$polling_place_id_2022), , drop = FALSE]
  frozen$match_id <- sprintf("AEC%05d", seq_len(nrow(frozen)))
  output <- data.frame(
    match_id = frozen$match_id, match_route = frozen$match_route, State = frozen$state_ab,
    PollingPlaceID_2019 = frozen$polling_place_id_2019,
    PollingPlaceID_2022 = frozen$polling_place_id_2022,
    PollingPlaceNm_2019 = frozen$polling_place_name_2019,
    PollingPlaceNm_2022 = frozen$polling_place_name_2022,
    NormalizedName_2019 = frozen$name_normalised_2019,
    NormalizedName_2022 = frozen$name_normalised_2022,
    Latitude_2019 = ifelse(frozen$coordinate_valid_2019, frozen$latitude_2019, NA_real_),
    Longitude_2019 = ifelse(frozen$coordinate_valid_2019, frozen$longitude_2019, NA_real_),
    Latitude_2022 = ifelse(frozen$coordinate_valid_2022, frozen$latitude_2022, NA_real_),
    Longitude_2022 = ifelse(frozen$coordinate_valid_2022, frozen$longitude_2022, NA_real_),
    MatchDistanceKm = frozen$match_distance_km,
    ExactIDOutlierAdjudication = frozen$exact_id_outlier_adjudication,
    X_2019 = frozen$X_2019, M_2019 = frozen$M_2019, Y_2019 = frozen$Y_2019,
    X_2022 = frozen$X_2022, M_2022 = frozen$M_2022, Y_2022 = frozen$Y_2022,
    stringsAsFactors = FALSE
  )
  if (file.exists(frozen_path) || file.exists(exclusion_path)) {
    stop("Refusing to overwrite an existing frozen match or exclusion file", call. = FALSE)
  }
  dir.create(dirname(frozen_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(output, frozen_path, row.names = FALSE, quote = TRUE)
  write.csv(exclusions, exclusion_path, row.names = FALSE, quote = TRUE)
  Sys.chmod(c(frozen_path, exclusion_path), mode = "0444", use_umask = FALSE)
  list(frozen = output, exclusions = exclusions,
       candidate_matches = nrow(candidates), exact_candidates = nrow(exact),
       secondary_candidates = nrow(secondary), denominator_exclusions = nrow(denominator_exclusions))
}
