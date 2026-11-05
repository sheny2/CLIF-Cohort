# Standalone CAPE COD-aligned early hydrocortisone analysis

This directory is self-contained. No script sources or executes a parent-folder
script. The existing CLIF-Cohort analyses are unchanged. Only the EHR source
outside this directory is required; all run outputs are inside `results/` here.

## Run on Discovery

From the existing CLIF-Cohort directory:

```bash
export TREATMENT_GRACE_HOURS=24
export CAPE_SEVERITY_MODE=pf
export CAPE_HYDRO_MODE=pragmatic_iv
export CLIF_DATA_CUTOFF_DATE=2026-10-31
export SUMMARY_SMALL_CELL_THRESHOLD=10
export CCW_MAX_PATIENTS=0
Rscript CAPE_COD_aligned/run_pipeline.R ../rclif_2_1
```

The source defaults to `../../rclif_2_1` relative to this new directory, equivalent
to `../rclif_2_1` from CLIF-Cohort. The runner works from any working directory
when no data argument is passed. R packages: arrow, data.table, dplyr, ggplot2.
A unique timestamped run directory is created by default. Optionally set
`CAPE_RUN_ID=primary_24h`; an existing run is refused. Each stage writes a log;
failure stops the pipeline. `COMPLETE.txt` is written only after every stage passes.
Use a new run ID for 48/96-hour sensitivity analyses. Twenty-four hours is primary;
changing grace changes the intervention, not just an estimator tuning parameter.

## Target and limits of emulation

The primary target is the 28-day mortality contrast between:

1. First qualifying IV hydrocortisone administration within 24 hours of the first
   simultaneous ICU/respiratory eligibility record.
2. No systemic steroid initiation during that same window.

After the grace window, usual care (including rescue steroids) is allowed in
both arms. Other systemic steroid regimens during grace are protocol deviations:
follow-up is censored at their first recorded administration, using the existing
natural-censoring weighting machinery. This comparator is deliberately an
**early initiation policy**, not the parent's never-treated-through-day-28 policy.
The treatment clone retains deaths during the final grace interval; surviving
noninitiators leave its risk set after that interval. Both clones enter at time
zero, and the main CCW analysis does not condition on surviving the grace period.

This is a pragmatic CAPE COD-aligned **initiation** analysis. It does NOT identify
the RCT's randomized assignment effect or enforce its full 4/7-day initial course,
8/14-day taper, or ICU-discharge stopping rule. No patient must survive or complete
future doses to be classified as an initiator. Before claiming exact protocol
emulation, review the RCT supplement's rescue, eligibility and regimen details.

RCT reference: Dequin et al., NEJM 2023, doi:10.1056/NEJMoa2215145,
https://www.nejm.org/doi/full/10.1056/NEJMoa2215145 .

## Cohort construction

- Combine overlapping hospitalizations and transfers with gaps <=6 hours using
  running maximum discharge times. Age >=18 and <120; one randomly selected
  eligible block per patient (seed 32284, configurable by COHORT_SELECTION_SEED).
- POA ICD-10 pneumonia J12-J18 in any diagnosis position; influenza, aspiration
  and VAP exclusions inherited from the clinical team's definition. Additional
  baseline codes exclude COVID (U07.1/J12.82), cystic fibrosis, tuberculosis and
  selected fungal/Pneumocystis infection codes. ICD coding is a proxy for clinical
  and radiological CAP, not proof of the RCT diagnostic criteria.
- First observed simultaneous ICU occupancy and qualifying respiratory record
  within 48 hours of block admission. The 48-hour support window is an explicit
  restriction to an early respiratory-failure subset, not a verbatim RCT criterion.
  Nonrebreather-mask and PSI-only pathways are omitted; CPAP-only records are not
  added to the source's IMV/NIPPV categories. Patients whose earlier eligibility
  was unrecorded can have imperfectly measured time zero.
- IMV/NIPPV requires recorded PEEP >=5 cm H2O within the preceding four hours.
  HFNO requires FiO2 >=0.50 and P/F <300 using arterial PO2 within the preceding
  four hours. No future measurement is used in these eligibility tests.
- Require ICU location at time zero, not merely a future ICU admission. No
  additional ED restriction is imposed, but the delivered source remains
  ED/ever-respiratory-support filtered; omitted admissions cannot be recovered.
- Exclude DNR/DNI-equivalent status, recorded prior death, and the clinical team's
  baseline shock proxy: any NE/epinephrine on HFNO/NIPPV; NE/epinephrine >0.25
  mcg/kg/min on IMV; active vasopressin; lactate >=4. Missing shock information is
  retained and flagged, so absence of recorded shock does not prove its absence.
- Exclude ANY observed positive systemic steroid administration from block
  admission through time zero, before one-block selection. This is an inpatient
  new-user proxy; outpatient/chronic exposure and competing clinical indications
  cannot be reliably excluded from this extract alone. Prior use outside the
  selected block is not assumed absent. An administration exactly at time zero
  is conservatively considered prebaseline exposure.
- Other RCT exclusions (e.g., postobstructive pneumonia, prolonged preadmission
  antibiotics, all competing indications) are not asserted to be recovered.

The four-hour lookback does not guarantee that independently sampled FiO2/PO2
represent identical physiological conditions. Validate source units and sample
timing before clinical interpretation. Death follow-up starts at EHR eligibility,
not the RCT randomization timestamp, which was unavailable in the EHR.

## Medication sources and definitions

`02_define_hydrocortisone.R` reads BOTH continuous and intermittent MAR tables.
Systemic steroid identities and IV/enteral routes are normalized; known held,
stopped, not-given, discontinued and ordered-only records are excluded. Positive
records without an action field are treated as administrations; if an action
field exists but is missing/unrecognized, positive-dose evidence is retained.
Verify the MAR action mapping with the source team. Continuous records without
a route column are assumed IV. Rate records are not integrated into daily totals.

Default `CAPE_HYDRO_MODE=pragmatic_iv` accepts:

- Continuous IV hydrocortisone at a recorded rate equivalent to 180-240 mg/day
  (a declared tolerance around 200 mg/day; units mg/hr, mg/day, mcg/min supported).
- Intermittent IV hydrocortisone 50-100 mg per administration. This supports common
  bolus practice but does NOT verify a sustained 200-mg/day course. A single bolus
  or a loading dose qualifies. Inspect the sensitive medication audit before use.

`CAPE_HYDRO_MODE=continuous_only` excludes intermittent doses from qualifying
initiation (they become incompatible regimens during grace). Positive systemic
steroid records with unsupported units or off-target doses are retained as
nonqualifying exposure and may cause deviation censoring, not untreated status.
Death/discharge/cutoff bounds are applied before choosing an initiation time.

## Schema options

Preflight writes schema and lab categories to the sensitive input-check folder.
Missing required tables/columns fail with an explicit message.

- `CAPE_PEEP_COLUMN=peep_set` can be changed to the actual PEEP field.
- `CAPE_PAO2_CATEGORIES='po2_arterial;pao2;pao2_arterial;arterial_po2'` must match
  arterial measurements only. P/F assumes CLIF standardized PO2 in mmHg.
- Default `CAPE_SEVERITY_MODE=pf` requires a configured arterial PO2 category.
- Explicit `CAPE_SEVERITY_MODE=sf_proxy` replaces HFNO P/F with S/F <315 and
  SpO2 80-97%. It remains a proxy, not equivalent RCT eligibility. PEEP is still
  required for IMV/NIPPV. Mode is recorded in run provenance and cohort design.

## Adjustment, estimation and diagnostics

The copied, independent analysis library runs federated CCW, its pooled
site-stratified identity control, baseline-only landmark IPW, IPCW without
cloning and unweighted per-protocol comparison. CCW is primary. Landmark is
conditional on survival/observation at grace end and is not the same estimand.
Landmark does not adjust for post-treatment laboratory values in this version.

Time-varying adjustment includes recent lactate and active NE, epinephrine and
vasopressin indicators (six-hour medication-record lookback). Ongoing infusions
without a recent rate record can be missed; these are recorded-exposure proxies,
not reconstructed infusion intervals. The nonrespiratory SOFA proxy retains its
component sum with an added count of missing components; this does not make the
proxy a validated complete SOFA score.

Other covariate processing follows the copied pipeline: recent-value lookbacks,
plausibility checks, 0.5/99.5% winsorization, pooled median single imputation and
missingness indicators. pH includes the arterial-category mapping. Variables
observed in <1% of person-intervals or constant after processing are dropped.
This is not multiple imputation or a guarantee against informative missingness.
Treatment denominator models are hospital-specific; the numerator is shared.
Weights are truncated at their 1st/99th percentiles. Review probability warnings,
overlap, interval-level balance and sensitivity to truncation. The reported SEs
condition on fitted weights and truncation; no nuisance-aware bootstrap is run.

## Outputs (inside this folder)

`results/<run-id>/` contains:

- `run_design.txt` and `source_versions.csv`: source/configuration and checksums.
- `logs/`: per-stage logs; may contain disclosive counts, keep on Discovery.
- `sensitive/input_checks/`: schema/category audit.
- `sensitive/cohort/`: patient/block files, crosswalks, exclusions and attrition.
- `sensitive/treatment/`: systemic steroid records, qualifying hydrocortisone,
  assignment, deviations, medication-source/dose/unit audit.
- `sensitive/analysis/`: prepared person-interval data and fitted analysis objects.
- `shareable/cohort_eda/`: suppressed EDA, six-hour death/treatment timing through
  day 28, 24/96-hour comparison counts and plots.
- `shareable/ccw_results/`: risk/RR/RD estimates, comparison plots, weight and
  balance diagnostics, model specification, warnings, and federated identity check.

Sensitive files and logs stay on Discovery. Aggregate files remain subject to
normal project disclosure review. The six-hour plots count qualifying IV
hydrocortisone here, not all systemic steroid drugs. EDA is descriptive; it does
not apply the CCW protocol-deviation censoring used in causal estimation.

## Verification

```bash
Rscript CAPE_COD_aligned/tests/test_end_to_end.R
Rscript CAPE_COD_aligned/tests/test_ccw_policy.R
Rscript CAPE_COD_aligned/tests/test_medication_rules.R
```

The end-to-end test creates synthetic CLIF tables in a temporary directory,
copies this standalone folder's R scripts there, and executes every stage with
two hospitals. It verifies baseline steroid exclusion, both medication sources,
and finite estimates. Synthetic runs do not establish clinical exchangeability
or validate accuracy of source diagnosis, medication or death ascertainment.
The clinical dataset is not present in the local workspace; run it on Discovery.
