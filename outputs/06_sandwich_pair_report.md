# NfL sandwich pair compatibility report

This report uses pair-aware epitope assignment across the two experimentally validated antibodies. Scores are deterministic proxy metrics intended for ranking and handoff to structure-modeling tools.

| Metric | Value |
|---|---:|
| Antibody 1 | 7-H11-D3-2-C7 |
| Antibody 1 epitope | C_boundary_368_377 (368-377) |
| Antibody 1 binding proxy | 74.77 |
| Antibody 2 | 15-C12-H6 |
| Antibody 2 epitope | Cys322_anchor_316_331 (316-331) |
| Antibody 2 binding proxy | 79.99 |
| Epitope overlap ratio | 0.0 |
| Linear epitope gap aa | 36 |
| Clash proxy score | 0.0 |
| Sandwich compatibility proxy | 84.52 |

Recommended orientation:

- Capture antibody: `7-H11-D3-2-C7`
- Detection antibody: `15-C12-H6`

Reason: capture uses the antibody with the higher sequence developability proxy; detection antibody should be checked for label-site interference.

Interpretation: low epitope overlap and a low clash proxy support taking this pair forward into trimer co-folding or docking. Final pair selection should use actual Fab1:NfL:Fab2 structures and interface metrics when available.
