import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
import prepare_antibody_template_inputs as prep  # noqa: E402


AA1_TO_AA3 = {value: key for key, value in prep.AA3_TO_AA1.items()}


def _atom_line(serial, atom, resname, chain, number, icode, x, y, z, element):
    return (
        f"ATOM  {serial:5d} {atom:>4s} {resname:>3s} {chain}{number:4d}{icode:1s}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}  "
    )


def _write_complete_pdb(path, chain_records, *, c_n=1.33, remarks=False):
    lines = ["REMARK PDBinfo-LABEL:    1 H1"] if remarks else []
    serial = 1
    offset = 0.0
    for chain, residues in chain_records:
        for number, icode, aa in residues:
            resname = AA1_TO_AA3[aa]
            atoms = prep.EXPECTED_HEAVY_ATOMS[aa]
            coordinates = {"N": (offset, 0, 0), "CA": (offset + 0.45, 0.8, 0),
                           "C": (offset + 1.33, 0, 0), "O": (offset + 1.33, -1.0, 0)}
            side = 0
            for atom in sorted(atoms):
                if atom not in coordinates:
                    side += 1
                    coordinates[atom] = (offset + 0.45, 0.8 + side * 0.4, side * 0.2)
                xyz = coordinates[atom]
                lines.append(_atom_line(serial, atom, resname, chain, number, icode,
                                        *xyz, atom[0]))
                serial += 1
            offset += 1.33 + c_n
        lines.append("TER")
        offset += 10.0
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


class AntibodyTemplateInputsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fasta_path = PROJECT_ROOT / "validation" / "experimentally_validated_antibodies.fasta"
        cls.evidence_path = PROJECT_ROOT / "input" / "antibody_templates" / "chothia_numbering_evidence.json"
        cls.h_csv = PROJECT_ROOT / "input" / "antibody_templates" / "nfl_H.csv"
        cls.l_csv = PROJECT_ROOT / "input" / "antibody_templates" / "nfl_KL.csv"
        cls.fasta = prep._read_fasta(cls.fasta_path)
        cls.tables = prep._read_anarci_tables(cls.h_csv, cls.l_csv)
        cls.numbering = prep._read_numbering_evidence(cls.evidence_path)
        cls.contracts = prep._validate_contracts(cls.fasta, cls.tables, cls.numbering)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.paired = {}
        self.broad = {}
        for template_id in prep.TEMPLATE_IDS:
            contract = self.contracts[template_id]
            records = []
            for pdb_chain, chain_name in (("H", "VH"), ("L", "VL")):
                labels = [x[0] for x in contract["anarci_by_chain"][chain_name]]
                sequence = contract["neutral_vh" if chain_name == "VH" else "neutral_vl"]
                records.append((pdb_chain, [
                    (int(label.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")),
                     label[len(label.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")):], aa)
                    for label, aa in zip(labels, sequence, strict=True)
                ]))
            path = self.root / f"{template_id}.pdb"
            _write_complete_pdb(path, records, remarks=True)
            self.paired[template_id] = path
            broad_path = self.root / f"{template_id}.broad_hlt.pdb"
            broad_remarks = "\n".join(
                f"REMARK PDBinfo-LABEL: {index:4d} {name}"
                for index, name in enumerate(prep.CDR_ORDER, 1)
            )
            atom_text = "\n".join(
                line for line in path.read_text().splitlines()
                if not line.startswith("REMARK PDBinfo-LABEL:")
            )
            broad_path.write_text(broad_remarks + "\n" + atom_text + "\n")
            self.broad[template_id] = broad_path

        generic_vh = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
        generic_vl = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
        sequence = generic_vh + prep.GENERIC_LINKER_OBSERVED + generic_vl
        self.generic = self.root / "generic.pdb"
        _write_complete_pdb(self.generic, [("A", [(i, "", aa) for i, aa in enumerate(sequence, 1)])])

    def test_builds_exact_hlt_and_neutral_scfv_blocked_bundle(self):
        output = self.root / "bundle"
        summary = prep.prepare_antibody_template_inputs(
            paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
            official_broad_hlt_pdbs=self.broad,
            output_dir=output,
        )
        self.assertEqual(summary["execution_state"], "blocked_pending_human_review")
        self.assertFalse(summary["real_model_handoff_authorized"])
        evidence = json.loads((output / "antibody_template_evidence.json").read_text())
        self.assertFalse(evidence["known_positive_cdr_sequences_used"])
        expected_lengths = {
            prep.TEMPLATE_IDS[0]: [7, 6, 3, 11, 7, 9],
            prep.TEMPLATE_IDS[1]: [8, 5, 11, 11, 7, 9],
        }
        for template_id in prep.TEMPLATE_IDS:
            hlt = output / f"{template_id}.exact_chothia.hlt.pdb"
            remarks = [line for line in hlt.read_text().splitlines()
                       if line.startswith("REMARK PDBinfo-LABEL:")]
            observed = [sum(line.endswith(name) for line in remarks) for name in prep.CDR_ORDER]
            self.assertEqual(observed, expected_lengths[template_id])
            self.assertEqual(set(prep._parse_pdb(hlt)), {"H", "L"})
            scfv = prep._parse_pdb(output / f"{template_id}.neutral_seed.scfv.pdb", allowed_chains={"A"})["A"]
            sequence = "".join(r.amino_acid for r in scfv)
            contract = self.contracts[template_id]
            self.assertEqual(sequence, contract["neutral_vh"] + prep.LINKER + contract["neutral_vl"])
            self.assertEqual(scfv[len(contract["neutral_vh"])].amino_acid, "A")
            self.assertIn("CB", scfv[len(contract["neutral_vh"])].by_name)
            record = evidence["templates"][template_id]
            self.assertTrue(record["official_converter_broad_labels_discarded_and_normalized"])
            self.assertEqual(record["scfv_metrics"]["generic_linker_observed"], prep.GENERIC_LINKER_OBSERVED)
            self.assertFalse(record["scfv_metrics"]["linker_is_standard_GGGGS3"])
            self.assertEqual(record["scfv_metrics"]["generic_linker_coordinate_or_identity_changes"], [])
        fragment = json.loads((output / "antibody_template_manifest_fragment.blocked.json").read_text())
        self.assertEqual(fragment["execution_state"], "blocked_pending_human_review")
        self.assertFalse(fragment["real_model_handoff_authorized"])

    def test_rejects_known_positive_or_non_neutral_sequence(self):
        template_id = prep.TEMPLATE_IDS[0]
        source_id = template_id.removeprefix("template_")
        records = []
        for pdb_chain, chain_name in (("H", "VH"), ("L", "VL")):
            labels = [x[0] for x in self.contracts[template_id]["anarci_by_chain"][chain_name]]
            sequence = self.fasta[(source_id, chain_name)]
            records.append((pdb_chain, [(int(label.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")),
                                        label[len(label.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")):], aa)
                                       for label, aa in zip(labels, sequence, strict=True)]))
        _write_complete_pdb(self.paired[template_id], records)
        with self.assertRaisesRegex(prep.AntibodyTemplatePreparationError, "not the locked neutral seed"):
            prep.prepare_antibody_template_inputs(
                paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
                output_dir=self.root / "rejected",
            )

    def test_rejects_wrong_chothia_insertion_code(self):
        path = self.paired[prep.TEMPLATE_IDS[1]]
        text = path.read_text()
        # Rewrite every atom of H31A as H31B without changing sequence.
        changed = "\n".join(
            line[:26] + "B" + line[27:] if line.startswith("ATOM") and line[21] == "H"
            and line[22:26].strip() == "31" and line[26] == "A" else line
            for line in text.splitlines()
        ) + "\n"
        path.write_text(changed)
        with self.assertRaisesRegex(prep.AntibodyTemplatePreparationError, "numbering/insertions"):
            prep.prepare_antibody_template_inputs(
                paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
                output_dir=self.root / "rejected",
            )

    def test_rejects_bad_peptide_geometry(self):
        path = self.paired[prep.TEMPLATE_IDS[0]]
        lines = path.read_text().splitlines()
        for index, line in enumerate(lines):
            if line.startswith("ATOM") and line[21] == "H" and line[22:26].strip() == "2" and line[12:16].strip() == "N":
                lines[index] = line[:30] + f"{99.0:8.3f}" + line[38:]
                break
        path.write_text("\n".join(lines) + "\n")
        with self.assertRaisesRegex(prep.AntibodyTemplatePreparationError, "Peptide C-N distance"):
            prep.prepare_antibody_template_inputs(
                paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
                output_dir=self.root / "rejected",
            )

    def test_pinned_anarci_hash_and_overwrite_are_enforced(self):
        corrupt = self.root / "nfl_H.csv"
        corrupt.write_bytes(self.h_csv.read_bytes() + b"\n")
        with self.assertRaisesRegex(prep.AntibodyTemplatePreparationError, "hash mismatch"):
            prep.prepare_antibody_template_inputs(
                paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
                output_dir=self.root / "hash_rejected", heavy_anarci_csv=corrupt,
            )
        output = self.root / "bundle"
        prep.prepare_antibody_template_inputs(
            paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
            official_broad_hlt_pdbs=self.broad, output_dir=output,
        )
        before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()}
        with self.assertRaisesRegex(prep.AntibodyTemplatePreparationError, "Refusing to overwrite"):
            prep.prepare_antibody_template_inputs(
                paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
                official_broad_hlt_pdbs=self.broad, output_dir=output,
            )
        prep.prepare_antibody_template_inputs(
            paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
            official_broad_hlt_pdbs=self.broad, output_dir=output, overwrite=True,
        )
        after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()}
        self.assertEqual(before, after)

    def test_accepts_two_strictly_validated_prepared_scfvs_and_rejects_partial(self):
        prepared = {}
        for template_id in prep.TEMPLATE_IDS:
            contract = self.contracts[template_id]
            sequence = contract["neutral_vh"] + prep.LINKER + contract["neutral_vl"]
            path = self.root / f"{template_id}.prepared.scfv.pdb"
            _write_complete_pdb(
                path, [("A", [(index, "", aa) for index, aa in enumerate(sequence, 1)])]
            )
            prepared[template_id] = path
        output = self.root / "prepared_bundle"
        summary = prep.prepare_antibody_template_inputs(
            paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
            prepared_scfv_pdbs=prepared, output_dir=output,
        )
        evidence = json.loads((output / "antibody_template_evidence.json").read_text())
        for template_id in prep.TEMPLATE_IDS:
            record = evidence["templates"][template_id]
            self.assertEqual(record["scfv_coordinate_source"], "externally_refined_domain_graft")
            self.assertEqual(record["prepared_scfv_input"]["sha256"], hashlib.sha256(prepared[template_id].read_bytes()).hexdigest())
        with self.assertRaisesRegex(prep.AntibodyTemplatePreparationError, "partial input is forbidden"):
            prep.prepare_antibody_template_inputs(
                paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
                prepared_scfv_pdbs={prep.TEMPLATE_IDS[0]: prepared[prep.TEMPLATE_IDS[0]]},
                output_dir=self.root / "partial",
            )

    def test_real_neutral_anarci_headers_and_upstream_hash_binding(self):
        neutral_h = Path("/tmp/nfl_template_raw/neutral_chothia_H.csv")
        neutral_l = Path("/tmp/nfl_template_raw/neutral_chothia_KL.csv")
        if not neutral_h.is_file() or not neutral_l.is_file():
            self.skipTest("Real neutral ANARCI regression fixtures are unavailable")
        parsed = prep._read_neutral_anarci_tables(neutral_h, neutral_l)
        self.assertEqual(
            set(parsed),
            {(template_id, chain) for template_id in prep.TEMPLATE_IDS for chain in ("VH", "VL")},
        )
        provenance = self.root / "provenance.json"
        payload = {
            "schema": prep.UPSTREAM_PROVENANCE_SCHEMA,
            "input_artifact_sha256": {
                "paired_fv_pdbs": {key: hashlib.sha256(path.read_bytes()).hexdigest() for key, path in self.paired.items()},
                "prepared_scfv_pdbs": {},
                "generic_scfv_pdb": hashlib.sha256(self.generic.read_bytes()).hexdigest(),
                "official_broad_hlt_pdbs": {key: hashlib.sha256(path.read_bytes()).hexdigest() for key, path in self.broad.items()},
            },
        }
        provenance.write_text(json.dumps(payload))
        output = self.root / "bound"
        prep.prepare_antibody_template_inputs(
            paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
            official_broad_hlt_pdbs=self.broad,
            neutral_heavy_anarci_csv=neutral_h,
            neutral_light_anarci_csv=neutral_l,
            upstream_provenance_path=provenance, output_dir=output,
        )
        evidence = json.loads((output / "antibody_template_evidence.json").read_text())
        self.assertFalse(evidence["preparation_cli_model_or_network_execution_performed"])
        self.assertTrue(evidence["upstream_coordinate_generation_provenance"]["input_hashes_verified"])
        payload["input_artifact_sha256"]["generic_scfv_pdb"] = "0" * 64
        provenance.write_text(json.dumps(payload))
        with self.assertRaisesRegex(prep.AntibodyTemplatePreparationError, "does not exactly bind"):
            prep.prepare_antibody_template_inputs(
                paired_fv_pdbs=self.paired, generic_scfv_pdb=self.generic,
                official_broad_hlt_pdbs=self.broad,
                upstream_provenance_path=provenance, output_dir=self.root / "bad_bound",
            )


if __name__ == "__main__":
    unittest.main()
