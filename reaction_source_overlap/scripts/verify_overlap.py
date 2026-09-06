import json
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/raw/USPTO_MIT/MIT_separated"
OUT = ROOT / "artifacts/experiment/day1_go_stop"


def raw_parts(value):
    if pd.isna(value):
        return ()
    return tuple(sorted(x.strip() for x in str(value).split(".") if x.strip()))


def canonical_parts(value):
    values = []
    for text in raw_parts(value):
        mol = Chem.MolFromSmiles(text)
        values.append("!RAW:" + text if mol is None else Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    return tuple(sorted(values))


def record_key(row, normalizer):
    return (normalizer(row.REACTANT), normalizer(row.REAGENT), normalizer(row.PRODUCT))


def indexed_keys(frame, normalizer):
    result = {}
    for index, row in enumerate(frame.itertuples(index=False)):
        result.setdefault(record_key(row, normalizer), index)
    return result


def main():
    train = pd.read_csv(DATA / "train.csv")
    val = pd.read_csv(DATA / "val.csv")
    test = pd.read_csv(DATA / "test.csv")
    raw_train, raw_val = indexed_keys(train, raw_parts), indexed_keys(val, raw_parts)
    can_train, can_val = indexed_keys(train, canonical_parts), indexed_keys(val, canonical_parts)
    counts = {"raw_train": 0, "raw_val": 0, "canonical_train": 0, "canonical_val": 0, "canonical_both": 0}
    canonical_seen, examples = set(), []
    duplicates = 0
    for index, row in enumerate(test.itertuples(index=False)):
        raw = record_key(row, raw_parts)
        canonical = record_key(row, canonical_parts)
        counts["raw_train"] += raw in raw_train
        counts["raw_val"] += raw in raw_val
        in_train, in_val = canonical in can_train, canonical in can_val
        counts["canonical_train"] += in_train
        counts["canonical_val"] += in_val
        counts["canonical_both"] += in_train and in_val
        duplicates += canonical in canonical_seen
        canonical_seen.add(canonical)
        if (in_train or in_val) and len(examples) < 20:
            examples.append(
                {
                    "test_index": index,
                    "matched_split": "both" if in_train and in_val else "train" if in_train else "val",
                    "matched_index": can_train.get(canonical, can_val.get(canonical)),
                    "reactant": row.REACTANT,
                    "reagent": "" if pd.isna(row.REAGENT) else row.REAGENT,
                    "product": row.PRODUCT,
                }
            )
    result = {
        **counts,
        "canonical_union": counts["canonical_train"] + counts["canonical_val"] - counts["canonical_both"],
        "canonical_test_duplicates": duplicates,
    }
    (OUT / "overlap_verification.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    pd.DataFrame(examples).to_csv(OUT / "overlap_examples.csv", index=False)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
