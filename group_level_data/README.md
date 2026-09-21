# Group-level derived data

Aggregate quantities from the stage-04 model fits, sufficient to inspect and verify the
main results without access to the underlying subject-level data.

Nothing here is subject-level. Under the NIMH Data Archive Data Use Certification
governing the source data, subject-level derived quantities — including the per-subject
generative effective connectivity matrices — may not be redistributed. Qualified
investigators who obtain their own executed Data Use Certification can regenerate them
from the source data with the code in this repository.

| File | Contents |
|---|---|
| `Ceffgroup_low_anhedonia.csv` | Group-level GEC matrix, low-anhedonia group (232 × 232) |
| `Ceffgroup_high_anhedonia.csv` | Group-level GEC matrix, high-anhedonia group (232 × 232) |
| `regional_trophic_levels_group_means.csv` | Per-region trophic level, group mean and SD for each group plus the high-minus-low difference (232 rows) |
| `group_summary_measures.csv` | Group mean, SD and n for trophic coherence and the two model-fit metrics |

## Conventions

The GEC matrices are asymmetric, with `Ceff[i, j]` the directed influence **from** region
`j` **to** region `i` (row = target, column = source).

Region indices follow the combined parcellation used throughout, 1–32 subcortex
(Tian Scale-II) and 33–232 cortex (Schaefer-200, 7 networks). Note these are 1-based here;
the Python notebooks use 0-based indices 0–31 and 32–231.

Group sizes are 139 low-anhedonia and 71 high-anhedonia, the participants with a cluster
assignment and complete resting-state imaging.
