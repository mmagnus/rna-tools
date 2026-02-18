#!/usr/bin/env python3
"""Extract an mmCIF structure from an OpenFold3 JSON response."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import warnings


def parse_args() -> argparse.Namespace:
      parser = argparse.ArgumentParser(
          description="Select an mmCIF structure from an OpenFold3 inference result."
      )
      parser.add_argument("json_path", type=Path, help="Path to openfold3.json")
      parser.add_argument(
          "-o",
          "--output",
          type=Path,
          help=(
              "Where to write the mmCIF. Defaults to <json_path>.mmcif. "
              "Use '-' to stream to stdout."
          ),
      )
      parser.add_argument(
          "--input-id",
          help="Which input entry to use when the JSON contains multiple items.",
      )
      parser.add_argument(
          "--structure-name",
          help="Pick a specific structure by its 'name' field.",
      )
      parser.add_argument(
          "--rank-by",
          default="confidence_score",
          help=(
              "Score field used to choose a structure when --structure-name is not set. "
              "Use 'none' to skip ranking and take the first structure."
          ),
      )
      parser.add_argument(
          "--no-pdb",
          action="store_true",
          help="Skip writing a PDB copy (conversion requires Biopython).",
      )
      parser.add_argument(
          "--pdb-output",
          type=Path,
          help=(
              "Where to write the PDB file. Defaults to the CIF path with a .pdb suffix."
          ),
      )
      return parser.parse_args()


def load_json(path: Path) -> dict:
      try:
          text = path.read_text(encoding="utf-8")
      except FileNotFoundError:
          sys.exit(f"Could not find '{path}'.")
      try:
          return json.loads(text)
      except json.JSONDecodeError as exc:
          sys.exit(f"Failed to parse '{path}': {exc}")


def select_output(outputs, input_id: str | None) -> dict:
      if not outputs:
          sys.exit("No outputs found in the JSON payload.")
      if input_id is None:
          if len(outputs) == 1:
              return outputs[0]
          available = ", ".join(o.get("input_id", "<missing>") for o in outputs)
          sys.exit(
              "Multiple input entries found. Please provide --input-id. "
              f"Available ids: {available}"
          )
      for entry in outputs:
          if entry.get("input_id") == input_id:
              return entry
      sys.exit(f"No output entry with input_id '{input_id}' found.")


def select_structure(structures, name: str | None, rank_by: str | None) -> dict:
      if not structures:
          sys.exit("No structures_with_scores entries were present.")
      if name:
          for struct in structures:
              if struct.get("name") == name:
                  return struct
          available = ", ".join(s.get("name", "<unnamed>") for s in structures)
          sys.exit(f"Structure '{name}' not available. Options: {available}")
      rank_field = None if rank_by and rank_by.lower() == "none" else rank_by
      if rank_field:
          scored = []
          for struct in structures:
              score = struct.get(rank_field)
              try:
                  scored.append((float(score), struct))
              except (TypeError, ValueError):
                  continue
          if scored:
              scored.sort(key=lambda item: item[0], reverse=True)
              return scored[0][1]
      return structures[0]


def confirm_destination(destination: Path) -> None:
      """Ask before overwriting an existing file."""
      if not destination.exists():
          return
      answer = input(f"'{destination}' exists. Overwrite? [y/N] ").strip().lower()
      if answer not in {"y", "yes"}:
          sys.exit("Aborted: destination already exists.")


def convert_cif_to_pdb(cif_path: Path, pdb_path: Path) -> None:
      try:
          from Bio.PDB import MMCIFParser, PDBIO
          from Bio.PDB.MMCIF2Dict import MMCIF2Dict
          from Bio.PDB.PDBExceptions import PDBConstructionWarning
      except ImportError as exc:
          sys.exit(
              "Biopython is required for PDB conversion. Install it with 'pip install biopython'."
          )

      class OccupancyAwareParser(MMCIFParser):
          def get_structure(self, structure_id, filename):  # type: ignore[override]
              with warnings.catch_warnings():
                  if self.QUIET:
                      warnings.filterwarnings(
                          "ignore", category=PDBConstructionWarning
                      )
                  self._mmcif_dict = MMCIF2Dict(filename)
                  self._ensure_occupancy_defaults()
                  self._build_structure(structure_id)
                  self._structure_builder.set_header(self._get_header())
              return self._structure_builder.get_structure()

          def _ensure_occupancy_defaults(self) -> None:
              if "_atom_site.occupancy" in self._mmcif_dict:
                  return
              coords = self._mmcif_dict.get("_atom_site.Cartn_x")
              if coords is None:
                  raise KeyError("_atom_site.occupancy")
              # Populate occupancy with 1.0 so Biopython can proceed.
              self._mmcif_dict["_atom_site.occupancy"] = ["1.00"] * len(coords)

      parser = OccupancyAwareParser(QUIET=True)
      structure = parser.get_structure("structure", str(cif_path))
      confirm_destination(pdb_path)
      pdb_path.parent.mkdir(parents=True, exist_ok=True)
      pdb_io = PDBIO()
      pdb_io.set_structure(structure)
      pdb_io.save(str(pdb_path))
      print(f"Wrote PDB to {pdb_path}", file=sys.stderr)


def emit_structure(structure: dict, destination: Path | None) -> None:
      cif_text = structure.get("structure")
      if not isinstance(cif_text, str):
          sys.exit("Selected structure did not include mmCIF text.")
      if destination is None:
          sys.stdout.write(cif_text)
          if not cif_text.endswith("\n"):
              sys.stdout.write("\n")
          return
      confirm_destination(destination)
      destination.parent.mkdir(parents=True, exist_ok=True)
      destination.write_text(cif_text, encoding="utf-8")
      print(
          f"Wrote mmCIF for '{structure.get('name', '<unnamed>')}' to {destination}",
          file=sys.stderr,
      )


def main() -> None:
      args = parse_args()
      destination = args.output
      if destination is None:
          destination = args.json_path.with_name(args.json_path.stem + ".cif")
      elif str(destination) == "-":
          destination = None
      pdb_destination: Path | None = None
      if not args.no_pdb:
          if destination is None:
              sys.exit("PDB conversion requires writing the CIF to a file.")
          pdb_destination = args.pdb_output
          if pdb_destination is None:
              pdb_destination = destination.with_suffix(".pdb")
      payload = load_json(args.json_path)
      output_entry = select_output(payload.get("outputs", []), args.input_id)
      structure = select_structure(
          output_entry.get("structures_with_scores", []),
          args.structure_name,
          args.rank_by,
      )
      emit_structure(structure, destination)
      if pdb_destination is not None:
          convert_cif_to_pdb(destination, pdb_destination)

if __name__ == "__main__":
        main()
