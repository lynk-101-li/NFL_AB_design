# NfL Antigen Truncation Inference Rationale

This document records the reasoning chain behind the NfL truncation module in `NFL_AB_design`. It is intentionally written as a method narrative, not merely as a result summary.

## 1. Why truncation inference is needed

The antibody design workflow does not start from an arbitrary full-length NfL antigen. The project begins with an NfL-related species migrating at roughly 25-35 kDa under non-reducing conditions and roughly 6-12 kDa after DTT reduction. The preferred working model is an antiparallel tetramer assembled from two Cys322-Cys322 disulfide-linked dimers. A disulfide-linked NfL dimer plus an unknown binding partner remains a viable alternative.

This matters for antibody design. If the diagnostic assay detects a disease- or sample-processing-associated NfL fragment, then antibodies should be evaluated against the fragment region that is actually present and reproducible, not against poorly constrained full-length protein regions.

## 2. Biological premise

Canonical human NfL/NEFL has a single cysteine residue at Cys322. Each proposed disulfide-linked NfL dimer must therefore contain Cys322 in both monomeric fragments.

The reasoning is:

1. Non-reducing migration: approximately 25-35 kDa.
2. Migration after DTT: approximately 6-12 kDa.
3. The roughly fourfold shift supports, but does not prove, a tetrameric species.
4. The preferred model is `4M -> (M-S-S-M)2`: two disulfide-linked dimers associate in an antiparallel, staggered arrangement.
5. The dimer-dimer interface may combine a coiled-coil hydrophobic core, electrostatic interactions, and geometric/conformational complementarity.
6. A disulfide-linked dimer plus an unknown binding partner is not excluded.
7. A disulfide-linked NfL monomer fragment must include Cys322, which lies in rod/coil-2B.

## 3. Why cathepsin-like cleavage is considered

The upstream project hypothesis is that tissue protease activity, especially cathepsin-like proteolysis, could generate the observed NfL fragment. SnapGene can display annotated protein features, but it does not automatically predict cathepsin cleavage sites like it does for DNA restriction enzyme sites.

The code therefore implements a transparent cleavage-preference model:

- Cysteine cathepsin-like scoring, approximating Cathepsin B/L/S/K/V-like preferences.
- Aspartyl cathepsin-like scoring, approximating Cathepsin D/E-like preferences.

The goal is not to claim exact enzyme specificity. The goal is to assign reproducible boundary support to every peptide bond and then combine that support with mass and Cys322 constraints.

## 4. Scoring model

For every peptide bond in NfL, the workflow records the P4-P4' local context.

Cysteine cathepsin-like scoring emphasizes:

- hydrophobic or aromatic P2 residues;
- basic P1 residues;
- dibasic context where relevant;
- penalty for Pro at P1'.

Cathepsin D/E-like scoring emphasizes:

- hydrophobic or aromatic residues at P1 and P1';
- Leu/Phe/Tyr/Trp-rich cleavage boundaries;
- penalty for unfavorable P1' Pro.

The output is not a biochemical rate constant. It is a deterministic, auditable proxy score for boundary plausibility.

## 5. Fragment enumeration

After scoring all peptide bonds, the workflow enumerates candidate fragments by pairing plausible N-terminal and C-terminal boundaries.

A retained fragment must:

- include Cys322;
- have monomer mass between the configured lower and upper limits;
- have a theoretical monomer mass compatible with the broad 6-12 kDa reducing band;
- have boundary support at both ends;
- be compatible with the preferred NfL rod/coil-2B antigen region.

The default constraints are stored in:

```text
input/antigen_truncation/truncation_constraints.json
```

## 6. Why aa 280-377 becomes the working antigen region

The highest-priority region is supported by three independent constraints:

1. Mass: fragments around aa 280-377 produce an approximately 11 kDa monomer, compatible with the broad reducing-band range. Non-reducing SDS-PAGE migration is not treated as an exact theoretical oligomer-mass measurement.
2. Chemistry: all retained core fragments contain Cys322, the only cysteine that can explain a Cys322-Cys322 homodimer.
3. Boundary plausibility: boundaries around W279/F280, F280/K281, K281/S282, L375/L376, L376/N377 and N377/V378 have cathepsin-like support.

The core candidates are:

- NfL 280-375
- NfL 281-376
- NfL 282-377

These candidates define the antigen context for the antibody design workflow.

## 7. How truncation inference feeds antibody design

The antibody workflow uses the inferred fragment region to choose epitope windows:

- N-terminal boundary region near aa 279-282;
- Cys322 anchor region;
- C-terminal boundary region near aa 368-377;
- sliding windows across the primary fragment.

This prevents the antibody analysis from drifting into regions that are less relevant to the observed fragment.

## 8. What remains experimental

The truncation module is a computational inference layer. It should be validated by experiments such as:

- LC-MS/MS of the non-reducing 25-35 kDa and reducing 6-12 kDa bands;
- no-enzyme or semi-tryptic searches for neo-termini;
- detection of Cys322-containing peptides;
- disulfide-linked Cys322-Cys322 peptide evidence;
- C322S comparison under reducing and non-reducing SDS-PAGE;
- SEC-MALS or native mass spectrometry to distinguish tetramer stoichiometry from a dimer-plus-partner complex;
- cross-linking mass spectrometry if the noncovalent dimer-dimer interface remains ambiguous.

DTT establishes that disulfide chemistry contributes to complex stability, but DTT alone cannot prove tetramer stoichiometry or show that the dimer-dimer interface is purely hydrophobic.
- in vitro cathepsin digestion with inhibitor controls;
- immediate alkylation controls to rule out artificial oxidation.

Only after these experiments should the inferred fragment boundaries be treated as established biochemical facts.
