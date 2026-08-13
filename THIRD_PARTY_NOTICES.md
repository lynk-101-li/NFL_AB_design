# Third-party notices

The repository-level MIT license applies only to original project code and documentation unless a file or dependency states otherwise.

## AlphaFold Protein Structure Database

`input/structures/AF-P07196-F1-model_v6.pdb` is the AlphaFold Protein Structure Database prediction for human NEFL/P07196, model version 6. The cropped target and structure evidence under `input/structures/` are derived from that prediction.

- Source: <https://alphafold.ebi.ac.uk/entry/P07196>
- Data license: [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)
- AlphaFold Data copyright: Google DeepMind / DeepMind Technologies Limited, as stated by AlphaFold DB
- Database attribution and current citation guidance: <https://alphafold.ebi.ac.uk/>

These coordinates are theoretical predictions and are not experimental structures.

## External model repositories and assets

The directories under `third_party/` are Git submodules pinned to upstream commits; they are not relicensed by this project. The machine-readable inventory is `config/model_components.json`.

- RFantibody and IgGM: MIT code licenses.
- Germinal: Apache-2.0 repository code; PyRosetta and optional predictors have separate terms.
- tFold: PolyForm Noncommercial 1.0.0 at the pinned revision.
- IgFold: JHU Academic Software License at the pinned revision.
- ImmuneBuilder: BSD-3-Clause code license.
- AlphaFold 3: Apache-2.0 code; model parameters, outputs and databases have separate terms.
- Chai-1: Apache-2.0 at the pinned revision.
- Boltz-2: MIT at the pinned revision.

This repository does not redistribute model checkpoints, AlphaFold databases or parameters, PyRosetta, or other restricted assets.

Students must review the upstream terms before downloading or executing them. In particular, PyRosetta requires a separate Rosetta license and must not be redistributed through this repository.

The generated neutral template coordinate files under `input/template_structures/` retain machine-readable tool/version and input-hash provenance in `input/template_structures/antibody_template_evidence.json` and `input/antibody_templates/template_coordinate_generation_provenance.json`.
