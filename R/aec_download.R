.aec_jurisdictions <- c("NSW", "VIC", "QLD", "WA", "SA", "TAS", "ACT", "NT")

.aec_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

.aec_single_csv_link <- function(list_items, label_pattern, href_pattern, context) {
  texts <- trimws(vapply(list_items, rvest::html_text2, character(1)))
  hit <- grepl(label_pattern, texts, ignore.case = TRUE, perl = TRUE)
  if (sum(hit) != 1L) {
    stop(sprintf("Expected exactly one official menu item for %s; found %d", context, sum(hit)), call. = FALSE)
  }
  links <- rvest::html_elements(list_items[[which(hit)]], "a[href]")
  hrefs <- rvest::html_attr(links, "href")
  csv <- grepl("[.]csv$", hrefs, ignore.case = TRUE) &
    grepl(href_pattern, basename(hrefs), ignore.case = TRUE, perl = TRUE)
  if (sum(csv) != 1L) {
    stop(sprintf("Expected exactly one validated CSV href for %s; found %d", context, sum(csv)), call. = FALSE)
  }
  hrefs[which(csv)]
}

discover_aec_links <- function(manifest_path = "config/aec_manifest.csv") {
  .aec_require_packages(c("curl", "xml2", "rvest"))
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  required <- c("year", "election_id", "house_menu_url")
  if (!all(required %in% names(manifest)) || nrow(manifest) != 2L) {
    stop("AEC election manifest must contain exactly two rows and the required fields", call. = FALSE)
  }
  rows <- list()
  z <- 0L
  for (i in seq_len(nrow(manifest))) {
    year <- as.integer(manifest$year[i])
    election_id <- as.character(manifest$election_id[i])
    menu_url <- manifest$house_menu_url[i]
    menu_response <- NULL
    for (attempt in seq_len(3L)) {
      menu_response <- tryCatch(curl::curl_fetch_memory(menu_url), error = identity)
      if (!inherits(menu_response, "error")) break
    }
    if (inherits(menu_response, "error") || menu_response$status_code != 200L) {
      stop("Could not retrieve official AEC menu: ", menu_url, call. = FALSE)
    }
    document <- xml2::read_html(rawToChar(menu_response$content))
    list_items <- rvest::html_elements(document, "li")
    year_rows <- list()
    for (jurisdiction in .aec_jurisdictions) {
      label <- paste0("^First preferences by candidate by polling place\\s*-\\s*", jurisdiction, "(?:\\s|$)")
      basename_pattern <- paste0("^HouseStateFirstPrefsByPollingPlaceDownload-",
                                 election_id, "-", jurisdiction, "[.]csv$")
      href <- .aec_single_csv_link(
        list_items, label, basename_pattern,
        sprintf("%d first preferences %s", year, jurisdiction)
      )
      year_rows[[length(year_rows) + 1L]] <- data.frame(
        election_year = year, election_id = election_id, jurisdiction = jurisdiction,
        file_role = "first_preferences_polling_place", official_menu_url = menu_url,
        resolved_file_url = xml2::url_absolute(href, menu_url),
        official_basename = basename(href), stringsAsFactors = FALSE
      )
    }
    metadata_href <- .aec_single_csv_link(
      list_items, "^Polling places(?:\\s|$)",
      paste0("^GeneralPollingPlacesDownload-", election_id, "[.]csv$"),
      sprintf("%d polling places metadata", year)
    )
    year_rows[[length(year_rows) + 1L]] <- data.frame(
      election_year = year, election_id = election_id,
      jurisdiction = "NATIONAL_METADATA", file_role = "polling_places_metadata",
      official_menu_url = menu_url,
      resolved_file_url = xml2::url_absolute(metadata_href, menu_url),
      official_basename = basename(metadata_href), stringsAsFactors = FALSE
    )
    year_rows <- do.call(rbind, year_rows)
    if (nrow(year_rows) != 9L ||
        !setequal(year_rows$jurisdiction[year_rows$file_role == "first_preferences_polling_place"],
                  .aec_jurisdictions) ||
        sum(year_rows$file_role == "polling_places_metadata") != 1L) {
      stop("Official link discovery did not yield exactly eight jurisdictions and one metadata file for ",
           year, call. = FALSE)
    }
    for (j in seq_len(nrow(year_rows))) {
      z <- z + 1L
      rows[[z]] <- year_rows[j, , drop = FALSE]
    }
  }
  discovered <- do.call(rbind, rows)
  rownames(discovered) <- NULL
  if (nrow(discovered) != 18L || anyDuplicated(discovered[c("election_year", "official_basename")])) {
    stop("Official discovery must yield exactly 18 unique files", call. = FALSE)
  }
  discovered
}

.aec_content_type <- function(response) {
  type <- response$type
  if (length(type) == 1L && !is.null(type) && nzchar(type)) return(type)
  headers <- rawToChar(response$headers)
  hit <- grep("^content-type:", strsplit(headers, "\\r?\\n")[[1L]],
              ignore.case = TRUE, value = TRUE)
  if (length(hit)) trimws(sub("^[^:]+:", "", hit[[length(hit)]])) else NA_character_
}

download_aec_files <- function(discovered, raw_root = "data/raw/aec") {
  .aec_require_packages(c("curl", "digest"))
  required <- c("election_year", "election_id", "jurisdiction", "file_role",
                "official_menu_url", "resolved_file_url", "official_basename")
  if (!all(required %in% names(discovered)) || nrow(discovered) != 18L) {
    stop("Discovered AEC inventory must contain exactly 18 validated rows", call. = FALSE)
  }
  manifest_rows <- vector("list", nrow(discovered))
  for (i in seq_len(nrow(discovered))) {
    year_dir <- file.path(raw_root, as.character(discovered$election_year[i]))
    dir.create(year_dir, recursive = TRUE, showWarnings = FALSE)
    target <- file.path(year_dir, discovered$official_basename[i])
    if (file.exists(target)) {
      stop("Refusing to overwrite pre-existing raw AEC file: ", target, call. = FALSE)
    }
    temporary <- paste0(target, ".download")
    if (file.exists(temporary)) stop("Stale temporary download exists: ", temporary, call. = FALSE)
    response <- tryCatch(
      curl::curl_fetch_disk(discovered$resolved_file_url[i], temporary),
      error = function(e) {
        if (file.exists(temporary)) unlink(temporary)
        stop("AEC download failed for ", discovered$resolved_file_url[i], ": ",
             conditionMessage(e), call. = FALSE)
      }
    )
    if (response$status_code != 200L) {
      if (file.exists(temporary)) unlink(temporary)
      stop("AEC download returned HTTP ", response$status_code, " for ",
           discovered$resolved_file_url[i], call. = FALSE)
    }
    if (!file.rename(temporary, target)) {
      if (file.exists(temporary)) unlink(temporary)
      stop("Could not atomically install downloaded AEC file: ", target, call. = FALSE)
    }
    # Raw bytes are hashed once, immediately after successful acquisition.
    raw_hash <- digest::digest(target, algo = "sha256", file = TRUE)
    Sys.chmod(target, mode = "0444", use_umask = FALSE)
    manifest_rows[[i]] <- data.frame(
      election_year = discovered$election_year[i], election_id = discovered$election_id[i],
      jurisdiction = discovered$jurisdiction[i], file_role = discovered$file_role[i],
      official_menu_url = discovered$official_menu_url[i],
      resolved_file_url = discovered$resolved_file_url[i],
      official_basename = discovered$official_basename[i],
      local_path = gsub("\\\\", "/", target),
      download_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      file_size_bytes = file.info(target)$size,
      sha256 = raw_hash,
      http_status = response$status_code, content_type = .aec_content_type(response),
      acquisition_status = "SUCCESS", stringsAsFactors = FALSE
    )
  }
  manifest <- do.call(rbind, manifest_rows)
  rownames(manifest) <- NULL
  if (nrow(manifest) != 18L || any(manifest$acquisition_status != "SUCCESS")) {
    stop("AEC acquisition did not produce exactly 18 successful rows", call. = FALSE)
  }
  manifest
}
