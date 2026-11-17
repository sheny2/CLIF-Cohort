# Site labels used by this standalone analysis. Real site IDs are retained;
# absent IDs and placeholder labels do not define an analysis site.
is_mapped_site <- function(x) {
  label <- tolower(trimws(as.character(x)))
  !is.na(label) & nzchar(label) &
    !label %in% c("na", "n/a", "unknown", "unmapped", "unspecified",
                  "hopkins_unspecified") &
    !grepl("(^|[_ -])(unknown|unmapped|unspecified)$", label)
}

assert_mapped_sites <- function(x) {
  if (any(!is_mapped_site(x))) {
    stop("Unmapped-site records found in saved inputs. Rerun scripts 01--07 ",
         "in a fresh RUN_ROOT with the updated site exclusion before EDA, ",
         "local meta-analysis, or NCO diagnostics.", call. = FALSE)
  }
  invisible(TRUE)
}
