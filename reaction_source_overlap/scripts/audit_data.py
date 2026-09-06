import hashlib
import json
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
CSV_DIR = ROOT / "data/raw/USPTO_MIT/MIT_separated"
MAPPED_DIR = ROOT / "data/raw/USPTO_MIT_mapped/data"
OUT = ROOT / "artifacts/experiment/day1_go_stop"
PROCESSED = ROOT / "data/processed"


def canonical_parts(value, clear_maps=False):
    if pd.isna(value) or not str(value).strip():
        return [], True
    parts, valid = [], True
    for text in str(value).split("."):
        text = text.strip()
        if not text:
            continue
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            parts.append("!RAW:" + text)
            valid = False
            continue
        if clear_maps:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
        parts.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    return sorted(parts), valid


def digest(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def normalize_split(name):
    frame = pd.read_csv(CSV_DIR / f"{name}.csv")
    rows = []
    for idx, row in enumerate(frame.itertuples(index=False)):
        reactant, reactant_ok = canonical_parts(row.REACTANT)
        reagent, _ = canonical_parts(row.REAGENT)
        product, product_ok = canonical_parts(row.PRODUCT)
        key = json.dumps([reactant, reagent, product], separators=(",", ":"))
        rows.append(
            {
                "split": name,
                "source_index": idx,
                "source_row": getattr(row, "_0", idx),
                "reactant": row.REACTANT,
                "reagent": "" if pd.isna(row.REAGENT) else row.REAGENT,
                "product": row.PRODUCT,
                "input": getattr(row, "input", ""),
                "reactant_canonical": ".".join(reactant),
                "product_canonical": ".".join(product),
                "reactant_key": digest(json.dumps(reactant, separators=(",", ":"))),
                "product_key": digest(json.dumps(product, separators=(",", ":"))),
                "reaction_hash": digest(key),
                "valid_reactant_product": reactant_ok and product_ok,
            }
        )
    return pd.DataFrame(rows)


def bond_code(mol, first, second):
    if mol is None or first is None or second is None:
        return "0"
    bond = mol.GetBondBetweenAtoms(first, second)
    return "0" if bond is None else str(bond.GetBondType())


def environment(mol, atom_index, radius):
    if mol is None or atom_index is None:
        return "_"
    bonds = list(Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_index))
    atoms = {atom_index}
    for bond_index in bonds:
        bond = mol.GetBondWithIdx(bond_index)
        atoms.update((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    center = mol.GetAtomWithIdx(atom_index)
    isotope = center.GetIsotope()
    center.SetIsotope(999)
    try:
        return Chem.MolFragmentToSmiles(
            mol,
            atomsToUse=sorted(atoms),
            bondsToUse=bonds,
            canonical=True,
            isomericSmiles=True,
        )
    finally:
        center.SetIsotope(isotope)


def mapped_record(line):
    reaction, center = line.rstrip("\n").rsplit(" ", 1)
    reactant_text, product_text = reaction.split(">>")
    reactant_mol = Chem.MolFromSmiles(reactant_text)
    product_mol = Chem.MolFromSmiles(product_text)
    if reactant_mol is None or product_mol is None:
        return None
    reactant_maps = {a.GetAtomMapNum(): a.GetIdx() for a in reactant_mol.GetAtoms() if a.GetAtomMapNum()}
    product_maps = {a.GetAtomMapNum(): a.GetIdx() for a in product_mol.GetAtoms() if a.GetAtomMapNum()}
    for atom in reactant_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    for atom in product_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    signatures = {0: [], 1: []}
    for edit in center.split(";"):
        try:
            first, second = (int(x) for x in edit.split("-"))
        except ValueError:
            return None
        ri, rj = reactant_maps.get(first), reactant_maps.get(second)
        pi, pj = product_maps.get(first), product_maps.get(second)
        if (ri is None and pi is None) or (rj is None and pj is None):
            return None
        old, new = bond_code(reactant_mol, ri, rj), bond_code(product_mol, pi, pj)
        for radius in signatures:
            endpoints = sorted(
                [
                    environment(reactant_mol, ri, radius) + ">" + environment(product_mol, pi, radius),
                    environment(reactant_mol, rj, radius) + ">" + environment(product_mol, pj, radius),
                ]
            )
            signatures[radius].append("|".join(endpoints + [old + ">" + new]))
    left_parts, left_ok = canonical_parts(reactant_text, clear_maps=True)
    product_parts, product_ok = canonical_parts(product_text, clear_maps=True)
    return {
        "template_r0": digest(";".join(sorted(signatures[0]))),
        "template_r1": digest(";".join(sorted(signatures[1]))),
        "mapped_left": left_parts,
        "mapped_product": product_parts,
        "mapped_valid": left_ok and product_ok,
        "center_changes": len(signatures[0]),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    train = normalize_split("train")
    val = normalize_split("val")
    test = normalize_split("test")

    template_r0, template_r1 = Counter(), Counter()
    template_failed_train = 0
    with (MAPPED_DIR / "train.txt").open(encoding="utf-8") as handle:
        for line in handle:
            record = mapped_record(line)
            if record is None:
                template_failed_train += 1
            else:
                template_r0[record["template_r0"]] += 1
                template_r1[record["template_r1"]] += 1

    test_templates = []
    with (MAPPED_DIR / "test.txt").open(encoding="utf-8") as handle:
        for row, line in zip(test.itertuples(index=False), handle):
            record = mapped_record(line)
            if record is not None:
                reactant_parts = row.reactant_canonical.split(".") if row.reactant_canonical else []
                reactant_subset = all(reactant_parts.count(x) <= record["mapped_left"].count(x) for x in set(reactant_parts))
                product_match = (row.product_canonical.split(".") if row.product_canonical else []) == record["mapped_product"]
                if not (record["mapped_valid"] and reactant_subset and product_match):
                    record = None
            test_templates.append(record)

    test["template_r0"] = [x["template_r0"] if x else "" for x in test_templates]
    test["template_r1"] = [x["template_r1"] if x else "" for x in test_templates]
    test["template_frequency_r0"] = test["template_r0"].map(template_r0).fillna(0).astype(int)
    test["template_frequency_r1"] = test["template_r1"].map(template_r1).fillna(0).astype(int)
    test["center_changes"] = [x["center_changes"] if x else pd.NA for x in test_templates]
    test["template_failed"] = test["template_r1"].eq("")

    reference_hashes = set(train.reaction_hash) | set(val.reaction_hash)
    test["exact_reference_overlap"] = test.reaction_hash.isin(reference_hashes)
    test["normalized_duplicate"] = test.duplicated("reaction_hash", keep="first")
    test["included"] = test.valid_reactant_product & ~test.exact_reference_overlap & ~test.normalized_duplicate
    test["exposure_level"] = "protected_exact"
    test.loc[test.exact_reference_overlap, "exposure_level"] = "unknown_exposure"

    train[["source_index", "reactant_canonical", "product_canonical", "reactant_key", "product_key", "reaction_hash", "valid_reactant_product"]].to_parquet(PROCESSED / "train_reference.parquet", index=False)
    test.to_parquet(PROCESSED / "test_pool.parquet", index=False)
    test[["source_index", "reaction_hash", "valid_reactant_product", "normalized_duplicate", "exact_reference_overlap", "exposure_level", "included"]].to_csv(OUT / "exposure_audit.csv", index=False)

    summary = {
        "rows": {"train": len(train), "val": len(val), "test": len(test)},
        "valid_reactant_product": {"train": int(train.valid_reactant_product.sum()), "val": int(val.valid_reactant_product.sum()), "test": int(test.valid_reactant_product.sum())},
        "normalized_duplicates": {"train": int(train.duplicated("reaction_hash").sum()), "val": int(val.duplicated("reaction_hash").sum()), "test": int(test.normalized_duplicate.sum())},
        "test_overlap_train_val": int(test.exact_reference_overlap.sum()),
        "template_failed": {"train_mapped": template_failed_train, "test_aligned": int(test.template_failed.sum())},
        "test_included": int(test.included.sum()),
    }
    (OUT / "data_audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
