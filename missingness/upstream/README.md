# upstream/

Verbatim copies of the functions that built the primary cohort, kept here so the missingness pipeline is self-contained on the
workbench (only the `missingness/` folder is uploaded) and so changes to them are tracked with the code that depends on them.

| file | copied from | used by |
|---|---|---|
| `cci_function.R` | `../function/cci_function.R` (`cci_calc_func`, `healthy_cci_func`) | `pipeline/compliance_cohort.R` (CCI of the sensitivity cohort, and the check against `pre_final_cci.txt`); tests |
| `ind_descriptive_func.R` | `../function/ind_descriptive_func.R` (`filter_pid_func`) | tests only: `R/compliance.R` must select exactly the people `filter_pid_func` selects |

Copied at commit `bd34333`. Do not edit these files here: edit the original in `../function/` and copy it again.
`tests/testthat/test-compliance.R` fails if a copy and its original differ. They are not sourced by `load_missingness()`
(they attach tidyverse packages); the scripts and tests source them explicitly.
