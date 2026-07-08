# Antibody Template Inputs

This directory is reserved for antibody starting templates or framework backgrounds.

The template FASTA is not the final experimental result. It is a neutral placeholder showing the expected format for VH/VL templates that could seed a prospective design or perturbation library.

Header format:

```text
>{template_id}|{VH_or_VL}|{description}
```

Use `X` characters only as placeholders in template CDR regions. Before running a real prospective design campaign, replace placeholders with generated or experimentally sourced CDR sequences.

The final experimentally validated antibodies are stored separately under:

```text
NFL_AB_design/validation/experimentally_validated_antibodies.fasta
```

They are used as a validation/replay set, not as primary input templates.
