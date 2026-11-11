# Explore the complete delivered rclif_2_1 extract

Run from CLIF-Cohort on Discovery:

```bash
Rscript Dataset_EDA/eda_rclif_dataset.R ../rclif_2_1
```

Or supply a NEW output directory:

```bash
Rscript Dataset_EDA/eda_rclif_dataset.R \
  /projects/emulate_icu/rclif_2_1 \
  Dataset_EDA/results/whole_extract_first_run
```

The default output is a new timestamped directory under Dataset_EDA/results/.
Existing output directories are refused. No original dataset or cohort scripts
are modified. Packages: arrow, data.table, dplyr, ggplot2.

This scans every immediate `clif_*.parquet` file or Parquet dataset directory.
No adult/pneumonia/treatment/ED/respiratory-support ANALYTIC filter is applied.
The source itself remains a selected extract; this script does not recover
admissions or historical care omitted upstream. Non-CLIF CSV/PDF shipping files
are not scanned. Run counts can be compared with shipped_flow_chart.csv separately.

## What to open first

1. START_HERE.md: guide to interpreting denominators.
2. sample_size_summary.csv: unique patients in patient/hospitalization tables,
   unique admissions and source row counts. Patients in the patient table can
   differ from patients represented in the hospitalization table.
3. table_inventory.csv: all table sizes, column counts, Parquet fragments, exact
   distinct patient/admission IDs where present, and linked patient coverage.
4. timestamp_ranges.csv: observed min/max for every recognized timestamp/date
   field, including admission, discharge, measurement, administration and death.
5. variable_dictionary.csv: every field name, Arrow type and inferred role.
6. column_missingness.csv: exact NA counts, blank-string counts, percentages and
   non-finite numeric values, across the WHOLE table (not our 4,497-patient cohort).
7. category_frequencies.csv: up to 50 most frequent category/code/unit values
   per field, with exact counts for eligible low-cardinality fields.
8. numeric_distributions.csv: finite-value n, mean, SD, min/max, quantiles.
   Labs, vitals, medication doses and numerical assessments are separated by
   measurement/drug category and unit where available. No unit conversion is
   applied; unspecified units are labeled accordingly.
9. measurement_groups.csv: which measurements/drugs and units are present.
10. hospital_length_of_stay.csv: nonnegative encounter LOS distribution and a
    count of negative LOS rows. Encounters are not merged into hospital blocks.
11. hospital_sizes_first_adt.csv: hospital-level encounter and linked patient
    sizes, attributed to the earliest valid recorded ADT hospital. Encounters
    without a usable hospital/time record are not included in this breakdown.
12. table_linkage_coverage.csv: unique encounter IDs absent from the supplied
    hospitalization roster and conflicting patient links. These are not claims
    about records outside the extract.
13. plots/: monthly admission/discharge/death record counts and sampled numeric
    histograms for the 12 largest numeric groups in each table.

Rows in event tables represent measurements/administrations, NOT independent
patients. Category percentages use table rows as denominator. ADT hospital
frequencies are ADT-row frequencies; use hospital_sizes_first_adt.csv for
encounter attribution. Blank death timestamps normally include living patients;
their frequency does not quantify unknown vital status. Observed date ranges do
not establish continuous coverage, data refresh dates, or mortality ascertainment.
Generic/arterial/venous laboratory category distinctions are preserved, not pooled
by clinical meaning. This is a source inventory, not a validated phenotype map.

## Exact statistics versus sampled distributions

- Table counts, distinct IDs, missingness and timestamp extrema/counts use all
  scanned rows and are exact for a completed table.
- Numeric n, mean, SD, min/max also use all finite values via streaming moments.
- Numeric quantiles and histograms use bounded uniform priority reservoirs,
  sampled across ALL batches. These are not first-row samples. Each group retains
  at most 2,000 numeric observations by default; the report states sample size
  and whether all finite values fit in that sample.
- Histograms show sampled observation counts, not estimated patient totals.
- If a category field exceeds the configured distinct-value limit, its entire
  frequency table is omitted and flagged. A first-batch/top-value approximation
  is never mislabeled as an exact frequency table.
- If the numeric-group limit is reached, inventory flags it. Profiles then cover
  only the retained groups; increase the limit for exhaustive group summaries.
- Patient/admission IDs are held transiently for exact counts and linkage, but
  their VALUES are never written. Identifier names appear only as schema metadata.
  No row-level samples or free-text contents are exported.

## Memory, runtime and errors

Arrow and data.table run single-threaded. Tables are read in record batches.
Numeric reservoirs, category dictionaries and number of numeric groups are capped.
Exact ID sets and the hospitalization/ADT mapping use memory proportional to the
number of distinct IDs: this part is NOT constant-memory. Default safety limit is
2 million distinct IDs per set; a table exceeding it is marked failed rather than
silently reporting approximate counts. Full-dataset scans can take substantial time.
All columns are scanned to characterize missingness, including nonexported text
fields, so very large text/nested fields can still increase batch memory use.

Optional environment settings:

```bash
export EDA_BATCH_ROWS=65536          # lower to 16384 if necessary
export EDA_RESERVOIR_N=2000          # increase for more precise quantiles
export EDA_MAX_CATEGORIES=1000       # field cardinality before frequencies omitted
export EDA_MAX_NUMERIC_GROUPS=2000   # per-table measurement/drug/unit groups
export EDA_MAX_UNIQUE_IDS=2000000    # raise only if HPC memory permits
export SUMMARY_SMALL_CELL_THRESHOLD=10
```

Outputs are checkpointed after each table. Recoverable R errors are recorded in
table_errors.csv and processing continues; PARTIAL.txt plus a nonzero exit means
the result must not be interpreted as a complete dataset profile. COMPLETE.txt
means all discovered tables and available derived summaries completed. Native
library segmentation faults cannot be caught by R; inspect the last logged table.

Low-frequency category values/counts, numeric moments/quantiles for small groups,
and small histogram cells are suppressed. Overall/group counts, dates and
missingness remain exact. ALL outputs are
sensitive working summaries: keep them on Discovery pending usual disclosure
review. Suppression here is not a guarantee that outputs are approved for export.

Timestamp extrema use typed Arrow timestamps/dates when available. String date
fields support ISO-like timestamps, explicit Z/numeric UTC offsets and date-only
values; timestamps without an offset are assumed UTC. Unparseable nonmissing
strings are counted. Verify source timezone conventions before temporal analysis.

## Test without clinical data

```bash
Rscript Dataset_EDA/tests/test_dataset_eda.R
```

The test constructs small synthetic CLIF tables, scans them with tiny batches,
checks exact sample sizes, linkage, missingness, dates, unit-stratified numeric
statistics, high-cardinality handling and aggregate-only export, and renders plots.
The actual rclif_2_1 source is not present in this local workspace.
