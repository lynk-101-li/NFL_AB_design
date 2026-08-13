# NfL Antigen Truncation Inference

## Objective

Infer NfL rod/coil-2B antigen fragments compatible with a roughly 6-12 kDa reducing band and a roughly 25-35 kDa non-reducing complex.

## Constraints

- Protein sequence: canonical human NEFL/NfL, 543 aa.
- Cys322 is the only cysteine in canonical human NfL.
- The preferred biochemical model is an antiparallel, staggered tetramer made from two Cys322-Cys322 disulfide-linked dimers.
- A disulfide-linked dimer plus an unknown binding partner remains an alternative explanation.
- DTT supports a disulfide contribution but does not alone establish stoichiometry or a purely hydrophobic interface.
- Candidate monomer fragments must include Cys322.
- Boundary support is estimated from cathepsin-like cleavage preferences.

## Cleavage Model

Cysteine cathepsin-like scoring emphasizes hydrophobic/aromatic P2 preference, basic P1 support, dibasic context, and P1' Pro penalty. Cathepsin D/E-like scoring emphasizes hydrophobic or aromatic P1/P1' boundaries. Scores are deterministic substrate-preference proxies and require experimental validation.

## Summary

- Total peptide bonds scored: 542
- Medium/high cleavage candidates: 87
- reducing-band/Cys322-compatible fragment candidates: 1339

## Top Cleavage Sites

|cut_after_position|bond|context_P4_to_P4prime|domain_at_cut|cysteine_cathepsin_like_score|cathepsin_D_E_like_score|overall_priority|
|---|---|---|---|---|---|---|
|122|L122\|L123|EAEL\|LVLR|rod:coil1A|1.9|5.5|high|
|138|L138\|Y139|FRAL\|YEQE|rod:coil1B|3.2|5.5|high|
|432|Y432\|L433|TSSY\|LMST|tail|1.2|5.5|high|
|123|L123\|V124|AELL\|VLRQ|rod:coil1A|3.2|5.2|high|
|9|Y9\|Y10|YEPY\|YSTS|head|1.7|5.1|medium|
|222|F222\|L223|EISF\|LKKV|rod:coil1B|1.9|5.1|medium|
|279|W279\|F280|AEEW\|FKSR|rod:linker2|1.6|5.1|medium|
|311|L311\|L312|SRRL\|LKAK|rod:coil2B|2.1|5.1|medium|
|368|Y368\|L369|MARY\|LKEY|rod:coil2B|3.5|5.1|medium|
|375|L375\|L376|YQDL\|LNVK|rod:coil2B|1.6|5.1|medium|
|392|L392\|L393|YRKL\|LEGE|rod:coil2B|2.6|5.1|medium|
|442|Y442\|Y443|FPSY\|YTSH|tail|1.7|5.1|medium|

## Top Fragment Candidates

|fragment|length_aa|monomer_avg_mass_kDa|disulfide_homodimer_avg_mass_kDa|two_dimer_tetramer_avg_mass_kDa|N_terminal_cut|C_terminal_cut|combined_boundary_score|reducing_band_mass_fit_score|
|---|---|---|---|---|---|---|---|---|
|280-375|96|11.041|22.081|44.161|W279\|F280|L375\|L376|10.2|100.0|
|282-377|96|10.993|21.984|43.969|K281\|S282|N377\|V378|8.4|100.0|
|281-376|96|11.007|22.013|44.025|F280\|K281|L376\|N377|7.3|100.0|
|282-375|94|10.766|21.53|43.06|K281\|S282|L375\|L376|9.9|100.0|
|280-376|97|11.154|22.307|44.614|W279\|F280|L376\|N377|8.9|100.0|
|280-368|89|10.151|20.301|40.601|W279\|F280|Y368\|L369|10.2|100.0|
|280-377|98|11.269|22.535|45.07|W279\|F280|N377\|V378|8.7|100.0|
|280-374|95|10.928|21.854|43.709|W279\|F280|D374\|L375|8.7|100.0|
|280-370|91|10.393|20.783|41.567|W279\|F280|K370\|E371|9.3|100.0|
|281-375|95|10.894|21.786|43.573|F280\|K281|L375\|L376|8.6|100.0|

## Interpretation

The highest-priority fragments cluster around NfL aa 280-377, consistent with a rod/coil-2B fragment containing Cys322. These inferred fragments are passed into the antibody epitope and ranking modules.
