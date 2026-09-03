detect_aec_header <- function(lines, required_fields) {
  if (!is.character(lines) || !length(lines)) stop("lines must be a non-empty character vector", call. = FALSE)
  if (!is.character(required_fields) || !length(required_fields)) stop("required_fields must be non-empty", call. = FALSE)
  hits <- vapply(lines, function(line) all(vapply(required_fields, grepl, logical(1), x = line, fixed = TRUE)), logical(1))
  if (sum(hits) != 1L) {
    stop(sprintf("Expected exactly one header row containing %s; found %d",
                 paste(required_fields, collapse = " and "), sum(hits)), call. = FALSE)
  }
  which(hits)
}

normalise_place_name <- function(x) {
  if (!requireNamespace("stringi", quietly = TRUE)) stop("Package 'stringi' is required", call. = FALSE)
  x <- stringi::stri_trans_nfkd(as.character(x))
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- toupper(x)
  x <- gsub("&", " AND ", x, fixed = TRUE)
  x <- gsub("[^A-Z0-9]", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

parse_aec_files <- function(...) {
  parse_aec_inventory(...)
}

.aec_name_key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", trimws(x)))

.aec_resolve_column <- function(data, aliases, concept, required = TRUE) {
  keys <- .aec_name_key(names(data))
  hits <- which(keys %in% .aec_name_key(aliases))
  if (!length(hits)) {
    if (required) stop("Missing required AEC field for ", concept, call. = FALSE)
    return(NA_character_)
  }
  if (length(hits) > 1L) {
    reference <- as.character(data[[hits[1L]]])
    agree <- vapply(hits[-1L], function(i) identical(as.character(data[[i]]), reference), logical(1))
    if (!all(agree)) stop("Conflicting aliases for AEC field ", concept, call. = FALSE)
  }
  names(data)[hits[1L]]
}

.aec_read_csv <- function(path, file_role) {
  lines <- readLines(path, n = 30L, warn = FALSE, encoding = "UTF-8")
  required <- if (identical(file_role, "first_preferences_polling_place")) {
    c("StateAb", "PollingPlaceID")
  } else {
    c("PollingPlaceID", "PollingPlaceTypeID")
  }
  header_row <- detect_aec_header(lines, required)
  data <- read.csv(
    path, skip = header_row - 1L, check.names = FALSE, colClasses = "character",
    na.strings = NULL, stringsAsFactors = FALSE, fill = TRUE
  )
  names(data) <- sub("^\\ufeff", "", names(data))
  list(data = data, header_row = header_row, pre_header = lines[seq_len(header_row - 1L)])
}

.aec_result_aliases <- list(
  state_ab = c("StateAb"), division_id = c("DivisionID"),
  division_name = c("DivisionNm", "DivisionName"), polling_place_id = c("PollingPlaceID"),
  polling_place_name = c("PollingPlace", "PollingPlaceNm", "PollingPlaceName"),
  candidate_id = c("CandidateID"), candidate_surname = c("Surname", "CandidateSurname"),
  candidate_given_name = c("GivenNm", "GivenName"),
  party_ab = c("PartyAb", "PartyAbbreviation"), party_nm = c("PartyNm", "PartyName"),
  ordinary_votes = c("OrdinaryVotes")
)

.aec_metadata_aliases <- list(
  state_ab = c("StateAb", "State"), division_id = c("DivisionID"),
  division_name = c("DivisionNm", "DivisionName"), polling_place_id = c("PollingPlaceID"),
  polling_place_name = c("PollingPlaceNm", "PollingPlace", "PollingPlaceName"),
  polling_place_type_id = c("PollingPlaceTypeID"),
  premises_name = c("PremisesNm"), address_1 = c("PremisesAddress1"),
  address_2 = c("PremisesAddress2"), address_3 = c("PremisesAddress3"),
  suburb = c("PremisesSuburb"), premises_state_ab = c("PremisesStateAb"),
  postcode = c("PremisesPostCode"), latitude = c("Latitude", "Lat"),
  longitude = c("Longitude", "Long", "Lon"),
  polling_place_type_label = c("PollingPlaceType", "PollingPlaceTypeNm", "PollingPlaceTypeName")
)

.aec_map_columns <- function(data, aliases, optional = character()) {
  vapply(names(aliases), function(concept) {
    .aec_resolve_column(data, aliases[[concept]], concept, required = !concept %in% optional)
  }, character(1))
}

.aec_trim <- function(x) trimws(as.character(x))

.aec_parse_votes <- function(x, context) {
  original <- .aec_trim(x)
  missing <- !nzchar(original)
  cleaned <- gsub(",", "", original, fixed = TRUE)
  value <- suppressWarnings(as.numeric(cleaned))
  invalid <- !missing & (!is.finite(value) | value < 0 | abs(value - round(value)) > 1e-12)
  if (any(invalid)) stop("Invalid OrdinaryVotes value in ", context, call. = FALSE)
  value[missing] <- NA_real_
  value
}

.aec_schema_rows <- function(manifest_row, parsed, mapping) {
  concepts <- rep("", length(names(parsed$data)))
  for (concept in names(mapping)) {
    if (!is.na(mapping[[concept]])) concepts[names(parsed$data) == mapping[[concept]]] <- concept
  }
  data.frame(
    election_year = unname(manifest_row$election_year),
    jurisdiction = unname(manifest_row$jurisdiction),
    file_role = unname(manifest_row$file_role),
    official_basename = unname(manifest_row$official_basename),
    header_row = unname(parsed$header_row),
    column_position = seq_along(names(parsed$data)),
    original_column_name = names(parsed$data),
    canonical_concept = concepts,
    stringsAsFactors = FALSE
  )
}

.aec_canonical_result <- function(data, mapping, year, jurisdiction) {
  out <- data.frame(
    election_year = as.integer(year),
    jurisdiction = jurisdiction,
    state_ab = .aec_trim(data[[mapping[["state_ab"]]]]),
    division_id = .aec_trim(data[[mapping[["division_id"]]]]),
    division_name = .aec_trim(data[[mapping[["division_name"]]]]),
    polling_place_id = .aec_trim(data[[mapping[["polling_place_id"]]]]),
    polling_place_name = .aec_trim(data[[mapping[["polling_place_name"]]]]),
    candidate_id = .aec_trim(data[[mapping[["candidate_id"]]]]),
    candidate_surname = .aec_trim(data[[mapping[["candidate_surname"]]]]),
    candidate_given_name = .aec_trim(data[[mapping[["candidate_given_name"]]]]),
    party_ab = .aec_trim(data[[mapping[["party_ab"]]]]),
    party_nm = .aec_trim(data[[mapping[["party_nm"]]]]),
    ordinary_votes = .aec_parse_votes(data[[mapping[["ordinary_votes"]]]],
                                      paste(year, jurisdiction)),
    stringsAsFactors = FALSE
  )
  out
}

.aec_canonical_metadata <- function(data, mapping, year) {
  value <- function(concept) {
    column <- mapping[[concept]]
    if (is.na(column)) rep("", nrow(data)) else .aec_trim(data[[column]])
  }
  data.frame(
    election_year = as.integer(year), state_ab = value("state_ab"),
    division_id = value("division_id"), division_name = value("division_name"),
    polling_place_id = value("polling_place_id"), polling_place_name = value("polling_place_name"),
    polling_place_type_id = value("polling_place_type_id"),
    polling_place_type_label = value("polling_place_type_label"),
    premises_name = value("premises_name"), address_1 = value("address_1"),
    address_2 = value("address_2"), address_3 = value("address_3"),
    suburb = value("suburb"), premises_state_ab = value("premises_state_ab"),
    postcode = value("postcode"), latitude_text = value("latitude"),
    longitude_text = value("longitude"), stringsAsFactors = FALSE
  )
}

parse_aec_inventory <- function(manifest_path = "results/manifests/aec_acquisition_manifest.csv") {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  if (nrow(manifest) != 18L || any(manifest$acquisition_status != "SUCCESS")) {
    stop("AEC acquisition manifest must contain exactly 18 successful rows", call. = FALSE)
  }
  results <- list()
  metadata <- list()
  schema <- list()
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    parsed <- .aec_read_csv(row$local_path, row$file_role)
    if (row$file_role == "first_preferences_polling_place") {
      mapping <- .aec_map_columns(parsed$data, .aec_result_aliases)
      key <- paste(row$election_year, row$jurisdiction, sep = "_")
      results[[key]] <- .aec_canonical_result(parsed$data, mapping,
                                               row$election_year, row$jurisdiction)
    } else {
      mapping <- .aec_map_columns(parsed$data, .aec_metadata_aliases,
                                  optional = "polling_place_type_label")
      metadata[[as.character(row$election_year)]] <-
        .aec_canonical_metadata(parsed$data, mapping, row$election_year)
    }
    schema[[i]] <- .aec_schema_rows(row, parsed, mapping)
  }
  list(manifest = manifest, results = results, metadata = metadata,
       schema_inventory = do.call(rbind, schema))
}

.aec_anomaly_row <- function(year, jurisdiction, role, type, count, details, severity = "AUDIT") {
  data.frame(election_year = as.integer(year), jurisdiction = jurisdiction,
             file_role = role, anomaly_type = type, count = as.integer(count),
             severity = severity, details = details, stringsAsFactors = FALSE)
}

audit_aec_inventory <- function(inventory) {
  anomaly <- list()
  party_rows <- list()
  denominator_rows <- list()
  type_rows <- list()
  metadata_rows <- list()
  low_denominator_rows <- list()
  z <- pz <- dz <- tz <- mz <- lz <- 0L

  for (key in names(inventory$results)) {
    x <- inventory$results[[key]]
    year <- unique(x$election_year)
    jurisdiction <- unique(x$jurisdiction)
    metadata <- inventory$metadata[[as.character(year)]]
    metadata_keys <- paste(metadata$state_ab, metadata$polling_place_id, sep = "|")
    result_keys <- paste(x$state_ab, x$polling_place_id, sep = "|")
    is_informal <- !nzchar(x$party_ab) & toupper(x$party_nm) == "INFORMAL"
    is_formal_candidate <- !is_informal
    group_key <- paste(x$state_ab, x$division_id, x$polling_place_id, sep = "|")
    checks <- list(
      duplicate_exact_rows = duplicated(x) | duplicated(x, fromLast = TRUE),
      duplicate_candidate_rows = duplicated(x[c("state_ab", "division_id", "polling_place_id", "candidate_id")]) |
        duplicated(x[c("state_ab", "division_id", "polling_place_id", "candidate_id")], fromLast = TRUE),
      duplicate_candidate_party_polling_place = duplicated(x[c("state_ab", "division_id", "polling_place_id", "candidate_id", "party_ab", "party_nm")]) |
        duplicated(x[c("state_ab", "division_id", "polling_place_id", "candidate_id", "party_ab", "party_nm")], fromLast = TRUE),
      missing_polling_place_id = !nzchar(x$polling_place_id),
      missing_party_ab_formal_candidate = is_formal_candidate & !nzchar(x$party_ab),
      missing_party_nm_formal_candidate = is_formal_candidate & !nzchar(x$party_nm),
      explicit_informal_rows_excluded_from_formal_denominator = is_informal,
      missing_ordinary_votes = is.na(x$ordinary_votes),
      negative_votes = !is.na(x$ordinary_votes) & x$ordinary_votes < 0,
      noninteger_votes = !is.na(x$ordinary_votes) & abs(x$ordinary_votes - round(x$ordinary_votes)) > 1e-12,
      state_file_mismatch = x$state_ab != jurisdiction,
      result_id_absent_from_metadata = !result_keys %in% metadata_keys
    )
    for (type in names(checks)) {
      z <- z + 1L
      anomaly[[z]] <- .aec_anomaly_row(
        year, jurisdiction, "first_preferences_polling_place", type,
        if (type == "result_id_absent_from_metadata") length(unique(result_keys[checks[[type]]])) else sum(checks[[type]]),
        "Count is zero when the structural check passes.",
        if (type == "explicit_informal_rows_excluded_from_formal_denominator") "EXPECTED" else
          if (any(checks[[type]])) "REVIEW" else "PASS"
      )
    }

    party <- aggregate(
      list(candidate_rows = x$ordinary_votes, total_ordinary_votes = x$ordinary_votes),
      list(election_year = x$election_year, jurisdiction = x$jurisdiction,
           party_ab = x$party_ab, party_nm = x$party_nm),
      FUN = function(v) sum(!is.na(v))
    )
    vote_totals <- aggregate(
      x$ordinary_votes,
      list(election_year = x$election_year, jurisdiction = x$jurisdiction,
           party_ab = x$party_ab, party_nm = x$party_nm),
      FUN = function(v) sum(v, na.rm = TRUE)
    )
    names(vote_totals)[5L] <- "total_ordinary_votes"
    party$total_ordinary_votes <- vote_totals$total_ordinary_votes
    party$is_alp_identifier <- toupper(party$party_ab) == "ALP"
    pz <- pz + 1L
    party_rows[[pz]] <- party

    split_rows <- split(seq_len(nrow(x)), group_key)
    formal_split_rows <- lapply(split_rows, function(ii) ii[is_formal_candidate[ii]])
    totals <- vapply(formal_split_rows, function(ii) sum(x$ordinary_votes[ii]), numeric(1))
    candidate_counts <- lengths(formal_split_rows)
    alp_counts <- vapply(formal_split_rows, function(ii) sum(toupper(x$party_ab[ii]) == "ALP"), integer(1))
    informal_counts <- vapply(split_rows, function(ii) sum(is_informal[ii]), integer(1))
    group_state_id <- vapply(split_rows, function(ii) result_keys[ii[1L]], character(1))
    metadata_type <- metadata$polling_place_type_id[match(group_state_id, metadata_keys)]
    ordinary_group <- metadata_type == "1"
    dz <- dz + 1L
    denominator_rows[[dz]] <- data.frame(
      election_year = year, jurisdiction = jurisdiction,
      polling_place_groups = length(split_rows), raw_rows = nrow(x),
      formal_candidate_rows = sum(is_formal_candidate), informal_rows_excluded = sum(is_informal),
      minimum_candidate_rows = min(candidate_counts), maximum_candidate_rows = max(candidate_counts),
      missing_vote_groups = sum(vapply(formal_split_rows, function(ii) any(is.na(x$ordinary_votes[ii])), logical(1))),
      noninteger_or_negative_groups = sum(!is.finite(totals) | totals < 0 | abs(totals - round(totals)) > 1e-12),
      nonpositive_formal_denominator_groups = sum(totals < 2),
      ordinary_polling_place_groups = sum(ordinary_group),
      nonpositive_ordinary_formal_denominator_groups = sum(totals[ordinary_group] < 2),
      groups_without_exactly_one_informal_row = sum(informal_counts != 1L),
      groups_without_exact_alp = sum(alp_counts == 0L), groups_with_multiple_alp_rows = sum(alp_counts > 1L),
      denominator_definition = "sum OrdinaryVotes over every candidate row within state-division-polling-place",
      stringsAsFactors = FALSE
    )
    low <- which(ordinary_group & totals < 2)
    if (length(low)) {
      lz <- lz + 1L
      first <- vapply(split_rows[low], `[`, integer(1), 1L)
      low_denominator_rows[[lz]] <- data.frame(
        election_year = year, jurisdiction = jurisdiction,
        division_id = x$division_id[first], division_name = x$division_name[first],
        polling_place_id = x$polling_place_id[first],
        polling_place_name = x$polling_place_name[first],
        formal_ordinary_votes_total = totals[low],
        disposition = "EXCLUDE_BEFORE_PHASE5B_ANALYSIS_M_BELOW_2",
        stringsAsFactors = FALSE
      )
    }
    z <- z + 1L
    anomaly[[z]] <- .aec_anomaly_row(
      year, jurisdiction, "denominator_structural_audit",
      "ordinary_formal_denominator_below_2", length(low),
      "Explicitly listed in aec_denominator_low_groups.csv; exclude before any Phase 5B analysis.",
      if (length(low)) "REVIEW" else "PASS"
    )
  }

  for (year_name in names(inventory$metadata)) {
    x <- inventory$metadata[[year_name]]
    year <- as.integer(year_name)
    lat <- suppressWarnings(as.numeric(x$latitude_text))
    lon <- suppressWarnings(as.numeric(x$longitude_text))
    coord_missing <- !nzchar(x$latitude_text) | !nzchar(x$longitude_text) | is.na(lat) | is.na(lon)
    coord_zero <- (!is.na(lat) & lat == 0) | (!is.na(lon) & lon == 0)
    coord_impossible <- (!is.na(lat) & (lat < -90 | lat > 90)) |
      (!is.na(lon) & (lon < -180 | lon > 180)) | coord_zero
    x$latitude <- lat
    x$longitude <- lon
    x$coordinate_valid <- !coord_missing & !coord_impossible
    inventory$metadata[[year_name]] <- x

    types <- sort(unique(x$polling_place_type_id))
    for (type in types) {
      rows <- x$polling_place_type_id == type
      names_for_type <- sort(unique(x$polling_place_name[rows]))
      labels <- sort(unique(x$polling_place_type_label[rows & nzchar(x$polling_place_type_label)]))
      tz <- tz + 1L
      type_rows[[tz]] <- data.frame(
        election_year = year, polling_place_type_id = type,
        supplied_type_labels = if (length(labels)) paste(labels, collapse = "; ") else "NOT_SUPPLIED",
        row_count = sum(rows),
        ppvc_name_count = sum(grepl("PPVC|PRE.?POLL", x$polling_place_name[rows], ignore.case = TRUE)),
        mobile_or_team_name_count = sum(grepl("MOBILE|TEAM", x$polling_place_name[rows], ignore.case = TRUE)),
        deterministic_examples = paste(head(names_for_type, 25L), collapse = "; "),
        stringsAsFactors = FALSE
      )
    }

    state_id <- paste(x$state_ab, x$polling_place_id, sep = "|")
    duplicated_id_rows <- duplicated(state_id) | duplicated(state_id, fromLast = TRUE)
    by_id <- split(seq_len(nrow(x)), state_id)
    multiple_division_ids <- sum(vapply(by_id, function(ii) length(unique(x$division_id[ii])) > 1L, logical(1)))
    inconsistent_ids <- sum(vapply(by_id, function(ii) {
      length(unique(x$polling_place_name[ii])) > 1L ||
        length(unique(paste(x$latitude_text[ii], x$longitude_text[ii], sep = "|"))) > 1L
    }, logical(1)))
    ordinary <- x$polling_place_type_id == "1"
    mz <- mz + 1L
    metadata_rows[[mz]] <- data.frame(
      election_year = year, metadata_rows = nrow(x), ordinary_rule = "PollingPlaceTypeID == 1",
      ordinary_rows = sum(ordinary), ordinary_unique_state_ids = length(unique(state_id[ordinary])),
      missing_coordinates_all = sum(coord_missing), impossible_coordinates_all = sum(coord_impossible),
      missing_coordinates_ordinary = sum(coord_missing & ordinary),
      impossible_coordinates_ordinary = sum(coord_impossible & ordinary),
      duplicated_state_id_rows = sum(duplicated_id_rows),
      duplicated_state_ids = sum(lengths(by_id) > 1L),
      ids_in_multiple_divisions = multiple_division_ids,
      ids_with_inconsistent_name_or_coordinates = inconsistent_ids,
      stringsAsFactors = FALSE
    )

    result_year <- do.call(rbind, inventory$results[startsWith(names(inventory$results), paste0(year, "_"))])
    result_ids <- unique(paste(result_year$state_ab, result_year$polling_place_id, sep = "|"))
    metadata_without_results <- !state_id %in% result_ids
    metadata_path <- inventory$manifest$local_path[
      inventory$manifest$election_year == year &
        inventory$manifest$file_role == "polling_places_metadata"
    ]
    raw_metadata_lines <- readLines(metadata_path, warn = FALSE, encoding = "UTF-8")
    malformed_quote_count <- sum(grepl(
      "(?<=,)\"[^\"]+\"[[:space:]]+[^,\"]", raw_metadata_lines, perl = TRUE
    ))
    meta_checks <- c(
      missing_polling_place_id = sum(!nzchar(x$polling_place_id)),
      duplicated_state_id_rows = sum(duplicated_id_rows),
      ids_in_multiple_divisions = multiple_division_ids,
      inconsistent_id_records = inconsistent_ids,
      missing_coordinates = sum(coord_missing),
      impossible_coordinates = sum(coord_impossible),
      metadata_ids_without_first_preference_rows = length(unique(state_id[metadata_without_results])),
      malformed_raw_csv_quote_fields = malformed_quote_count
    )
    for (type in names(meta_checks)) {
      z <- z + 1L
      anomaly[[z]] <- .aec_anomaly_row(
        year, "NATIONAL_METADATA", "polling_places_metadata", type,
        meta_checks[[type]], "Reported metadata integrity/coverage diagnostic.",
        if (meta_checks[[type]] > 0L) "REVIEW" else "PASS"
      )
    }
  }

  party_inventory <- do.call(rbind, party_rows)
  rownames(party_inventory) <- NULL
  list(
    inventory = inventory,
    schema_inventory = inventory$schema_inventory,
    party_inventory = party_inventory,
    polling_place_type_inventory = do.call(rbind, type_rows),
    metadata_summary = do.call(rbind, metadata_rows),
    denominator_audit = do.call(rbind, denominator_rows),
    low_denominator_groups = if (length(low_denominator_rows)) do.call(rbind, low_denominator_rows) else data.frame(),
    anomaly_inventory = do.call(rbind, anomaly)
  )
}
