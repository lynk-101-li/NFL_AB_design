# NfL cathepsin site map for SnapGene

Input sequence: canonical human NEFL/NfL UniProt P07196, 543 aa.

SnapGene can import protein sequences as FASTA/GenPept/UniProt/SwissProt/GCG, but its enzyme display is for restriction enzyme sites on nucleotide sequences. For cathepsin/protease cuts, use this annotated GenPept file as custom protein features:

- `nfl_cathepsin_annotated_for_snapgene.gp`

Open in SnapGene with File -> Open. In Map/Sequence view, show Features. The file marks:

- UniProt domains, including coil 2B at aa 281-396
- Cys322 as the disulfide anchor candidate
- top Cys322-containing fragment candidates for the reduction-sensitive oligomer hypothesis
- medium/high cathepsin-like cleavage candidates around coil 2B

Companion CSV files:

- `nfl_all_peptide_bonds_cathepsin_scores.csv`: all 542 peptide bonds with substrate-specificity-based cathepsin-like scores
- `nfl_medium_high_cathepsin_candidate_sites.csv`: filtered medium/high candidate cuts
- `nfl_reduction_sensitive_fragment_candidates.csv`: Cys322-containing monomer candidates compatible with the 6-12 kDa reducing band; the non-reducing band is not treated as an exact dimer-mass constraint

Most likely starting hypotheses from the scoring table:

- aa 280-375: monomer about 11.04 kDa; disulfide homodimer about 22.08 kDa; boundaries W279|F280 and L375|L376
- aa 281-376: monomer about 11.01 kDa; disulfide homodimer about 22.01 kDa; boundaries F280|K281 and L376|N377
- aa 282-376 / 282-377: monomer about 10.88-10.99 kDa; disulfide homodimer about 21.76-21.99 kDa; N boundary K281|S282

Important caveat: cathepsin predictions are specificity-based candidates. Validate with non-reducing/reducing LC-MS/MS and semi-/no-enzyme searches for neo-termini.
