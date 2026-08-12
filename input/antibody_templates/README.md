# Antibody Template Inputs

This directory is reserved for neutral example antibody starting templates or framework backgrounds.

The template FASTA is not the final experimental result. It is a neutral placeholder showing the expected format for VH/VL templates that could seed a prospective design or perturbation library.

Header format:

```text
>{template_id}|{VH_or_VL}|{description}
```

Use `X` characters only as placeholders in template CDR regions. Before running a real prospective design campaign, replace placeholders with generated or experimentally sourced CDR sequences.

The two experimentally validated antibodies are stored separately under:

```text
NFL_AB_design/validation/experimentally_validated_antibodies.fasta
```

In the current campaign their VH/VL framework residues provide two distinct
framework sources only after all six CDRs are masked. Their known CDR amino-acid
identities are not generation features. The complete validated sequences are
introduced only after prospective ranking as explicitly labeled retrospective
positive controls.
