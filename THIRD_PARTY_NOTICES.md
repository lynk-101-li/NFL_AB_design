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

RFantibody, IgGM, Germinal, ImmuneBuilder, PyRosetta, AlphaFold parameters, Chai, ColabDesign and their model weights are third-party works governed by their own licenses and terms. This repository does not redistribute RFantibody, IgGM, Germinal, PyRosetta or AlphaFold parameter bundles.

Students must review the upstream terms before downloading or executing them. In particular, PyRosetta requires a separate Rosetta license and must not be redistributed through this repository.

The generated neutral template coordinate files under `input/template_structures/` retain machine-readable tool/version and input-hash provenance in `input/template_structures/antibody_template_evidence.json` and `input/antibody_templates/template_coordinate_generation_provenance.json`.
