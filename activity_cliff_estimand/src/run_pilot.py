import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Levenshtein import distance
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error

from prepare_data import cliff_pairs
from prepare_splits import mmp_cores


def ecfp(smiles, bits=2048):
    mols = [Chem.MolFromSmiles(value) for value in smiles]
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=bits) for mol in mols]
    array = np.zeros((len(fps), bits), dtype=np.uint8)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, array[i])
    return mols, fps, array


def nearest_distance(query_fps, train_fps):
    return np.array([1.0 - max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)) for fp in query_fps])


def scaffold_fps(mols):
    output = []
    for mol in mols:
        try:
            core = MakeScaffoldGeneric(mol)
        except Exception:
            core = GetScaffoldForMol(mol)
        output.append(AllChem.GetMorganFingerprintAsBitVect(core, 2, nBits=1024))
    return output


def pair_table(test, train_distance):
    smiles = test["canonical_smiles"].tolist()
    mols, morgan, _ = ecfp(smiles, bits=1024)
    scaffolds = scaffold_fps(mols)
    exp_nm = test["exp_nm"].to_numpy(float)
    weights = np.array([Descriptors.MolWt(mol) for mol in mols])
    core_sets = [mmp_cores(value) for value in smiles]
    core_sizes = {
        core: Chem.MolFromSmiles(core).GetNumHeavyAtoms()
        for cores in core_sets
        for core in cores
    }
    rows = []
    for i in range(len(test) - 1):
        tani = DataStructs.BulkTanimotoSimilarity(morgan[i], morgan[i + 1 :])
        scaf = DataStructs.BulkTanimotoSimilarity(scaffolds[i], scaffolds[i + 1 :])
        for offset, (tani_value, scaf_value) in enumerate(zip(tani, scaf), start=1):
            j = i + offset
            lev = 1.0 - distance(smiles[i], smiles[j]) / max(len(smiles[i]), len(smiles[j]))
            if tani_value < 0.9 and scaf_value < 0.9 and lev < 0.9:
                continue
            fold = max(exp_nm[i], exp_nm[j]) / min(exp_nm[i], exp_nm[j])
            shared_cores = core_sets[i] & core_sets[j]
            edit_size = np.nan
            if shared_cores:
                edit_size = (
                    mols[i].GetNumHeavyAtoms()
                    + mols[j].GetNumHeavyAtoms()
                    - 2 * max(core_sizes[core] for core in shared_cores)
                )
            rows.append(
                {
                    "i": i,
                    "j": j,
                    "pair_type": "cliff" if fold > 10.0 else "noncliff",
                    "max_similarity": max(tani_value, scaf_value, lev),
                    "avg_mw": (weights[i] + weights[j]) / 2,
                    "mw_diff": abs(weights[i] - weights[j]),
                    "train_distance_mean": (train_distance[i] + train_distance[j]) / 2,
                    "mmp_edit_size": edit_size,
                }
            )
    return pd.DataFrame(rows)


def match_pairs(pairs):
    cliffs = pairs[pairs["pair_type"] == "cliff"].reset_index(drop=True)
    noncliffs = pairs[pairs["pair_type"] == "noncliff"].reset_index(drop=True)
    if cliffs.empty or noncliffs.empty:
        return pd.DataFrame()
    columns = ["max_similarity", "avg_mw", "mw_diff", "train_distance_mean", "mmp_edit_size"]
    combined = pd.concat([cliffs[columns], noncliffs[columns]], ignore_index=True)
    scale = combined.std(ddof=0).replace(0, 1).fillna(1)
    a = cliffs[columns].fillna(combined.median()).to_numpy() / scale.to_numpy()
    b = noncliffs[columns].fillna(combined.median()).to_numpy() / scale.to_numpy()
    row_ind, col_ind = linear_sum_assignment(((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2))
    rows = []
    for match_id, (row, col) in enumerate(zip(row_ind, col_ind)):
        rows.append(
            {
                "match_id": match_id,
                "cliff_i": int(cliffs.iloc[row]["i"]),
                "cliff_j": int(cliffs.iloc[row]["j"]),
                "noncliff_i": int(noncliffs.iloc[col]["i"]),
                "noncliff_j": int(noncliffs.iloc[col]["j"]),
            }
        )
    return pd.DataFrame(rows)


def rmse(y, prediction):
    return float(np.sqrt(mean_squared_error(y, prediction)))


def indexed_rmse(y, prediction, indices):
    indices = sorted(set(indices))
    return rmse(y[indices], prediction[indices]) if indices else np.nan


def risk_metrics(y, prediction, uncertainty):
    order = np.argsort(uncertainty, kind="stable")
    coverage = np.arange(1, len(y) + 1) / len(y)
    risks = np.array([rmse(y[order[:k]], prediction[order[:k]]) for k in range(1, len(y) + 1)])
    keep = max(1, int(np.floor(0.8 * len(y))))
    rho = spearmanr(uncertainty, np.abs(y - prediction)).statistic
    return coverage, risks, float(np.trapz(risks, coverage)), float(risks[keep - 1] / risks[-1]), float(rho)


def load_features(representation, dataset, frame, feature_dir, config):
    smiles = frame["canonical_smiles"].tolist()
    _, fps, matrix = ecfp(smiles, bits=config["prediction_ecfp"]["bits"])
    if representation == "ecfp":
        return matrix, fps
    cache = np.load(Path(feature_dir) / f"{dataset}.npz")
    lookup = {smiles_value: embedding for smiles_value, embedding in zip(cache["smiles"], cache["embeddings"])}
    return np.vstack([lookup[value] for value in smiles]), fps


def run_one(dataset, environment, frame, representation, feature_dir, config, matching_source=None):
    features, fps = load_features(representation, dataset, frame, feature_dir, config)
    y = frame["activity"].to_numpy(float)
    train_idx = np.where(frame["split"].to_numpy() == "train")[0]
    cal_idx = np.where(frame["split"].to_numpy() == "calibration")[0]
    test_idx = np.where(frame["split"].to_numpy() == "test")[0]
    mean, std = y[train_idx].mean(), y[train_idx].std(ddof=0)
    predictions = []
    cal_predictions = []
    prediction_rows = []
    for seed in config["seeds"]:
        model = RandomForestRegressor(random_state=seed, **config["rf"])
        model.fit(features[train_idx], (y[train_idx] - mean) / std)
        test_prediction = model.predict(features[test_idx]) * std + mean
        cal_prediction = model.predict(features[cal_idx]) * std + mean
        predictions.append(test_prediction)
        cal_predictions.append(cal_prediction)
        for split_name, indices, values in (("calibration", cal_idx, cal_prediction), ("test", test_idx, test_prediction)):
            for index, value in zip(indices, values):
                prediction_rows.append(
                    {
                        "dataset": dataset,
                        "environment": environment,
                        "representation": representation,
                        "seed": seed,
                        "split": split_name,
                        "canonical_smiles": frame.iloc[index]["canonical_smiles"],
                        "y_true": y[index],
                        "prediction": value,
                    }
                )
    predictions = np.vstack(predictions)
    cal_predictions = np.vstack(cal_predictions)
    ensemble_mean = predictions.mean(axis=0)
    ensemble_variance = predictions.var(axis=0, ddof=1)
    test_distance = nearest_distance([fps[i] for i in test_idx], [fps[i] for i in train_idx])
    y_test = y[test_idx]
    test = frame.iloc[test_idx].reset_index(drop=True)
    if matching_source is None:
        pairs = pair_table(test, test_distance)
        matches = match_pairs(pairs)
        cliff_index_pairs = pairs.loc[pairs["pair_type"] == "cliff", ["i", "j"]].to_numpy()
    else:
        matches = matching_source.reset_index(drop=True)
        cliff_index_pairs = np.asarray(cliff_pairs(test), dtype=int)
    cliff_indices = cliff_index_pairs.ravel().tolist()
    matched_cliff = matches[["cliff_i", "cliff_j"]].to_numpy().ravel().tolist() if not matches.empty else []
    matched_noncliff = matches[["noncliff_i", "noncliff_j"]].to_numpy().ravel().tolist() if not matches.empty else []
    metrics = {
        "dataset": dataset,
        "environment": environment,
        "representation": representation,
        "n_train": len(train_idx),
        "n_calibration": len(cal_idx),
        "n_test": len(test_idx),
        "rmse": rmse(y_test, ensemble_mean),
        "mae": float(mean_absolute_error(y_test, ensemble_mean)),
        "cliff_rmse": indexed_rmse(y_test, ensemble_mean, cliff_indices),
        "matched_cliff_rmse": indexed_rmse(y_test, ensemble_mean, matched_cliff),
        "matched_noncliff_rmse": indexed_rmse(y_test, ensemble_mean, matched_noncliff),
        "n_cliff_pairs": int(len(cliff_index_pairs)),
        "n_matched_pairs": len(matches),
    }
    metrics["cliff_ratio"] = metrics["matched_cliff_rmse"] / metrics["matched_noncliff_rmse"]
    risk_rows = []
    uq_rows = []
    for method, values in (("ensemble_variance", ensemble_variance), ("ecfp_distance", test_distance)):
        coverage, risks, aurc, sr80, rho = risk_metrics(y_test, ensemble_mean, values)
        metrics[f"{method}_aurc"] = aurc
        metrics[f"{method}_sr80"] = sr80
        metrics[f"{method}_spearman_abs_error"] = rho
        for item, uncertainty, error in zip(test["canonical_smiles"], values, np.abs(y_test - ensemble_mean)):
            uq_rows.append(
                {
                    "dataset": dataset,
                    "environment": environment,
                    "representation": representation,
                    "canonical_smiles": item,
                    "method": method,
                    "uncertainty": uncertainty,
                    "absolute_error": error,
                }
            )
        for value, risk in zip(coverage, risks):
            risk_rows.append(
                {
                    "dataset": dataset,
                    "environment": environment,
                    "representation": representation,
                    "method": method,
                    "coverage": value,
                    "rmse": risk,
                }
            )
    cal_scores = np.abs(y[cal_idx] - cal_predictions.mean(axis=0))
    rank = min(len(cal_scores), int(np.ceil((len(cal_scores) + 1) * (1 - config["conformal_alpha"]))))
    quantile = float(np.partition(cal_scores, rank - 1)[rank - 1])
    metrics["conformal_coverage_90"] = float(np.mean(np.abs(y_test - ensemble_mean) <= quantile))
    metrics["conformal_mean_width"] = 2 * quantile
    matches = matches.assign(dataset=dataset, environment=environment, representation=representation)
    return metrics, prediction_rows, uq_rows, risk_rows, matches


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--split-manifest", required=True)
    parser.add_argument("--representation", choices=["ecfp", "molformer"], required=True)
    parser.add_argument("--feature-dir")
    parser.add_argument("--matching-source")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    manifest = pd.read_csv(args.split_manifest)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    matching_source = pd.read_csv(args.matching_source) if args.matching_source else None
    metrics, predictions, uq_rows, risk_rows, matches = [], [], [], [], []
    for (dataset, environment), frame in manifest.groupby(["dataset", "environment"], sort=True):
        source = None
        if matching_source is not None:
            source = matching_source[
                (matching_source["dataset"] == dataset) & (matching_source["environment"] == environment)
            ].copy()
        result = run_one(
            dataset,
            environment,
            frame.reset_index(drop=True),
            args.representation,
            args.feature_dir,
            config,
            source,
        )
        metrics.append(result[0])
        predictions.extend(result[1])
        uq_rows.extend(result[2])
        risk_rows.extend(result[3])
        matches.append(result[4])
        print(f"DONE {dataset} {environment} {args.representation} RMSE={result[0]['rmse']:.4f}", flush=True)
        pd.DataFrame(metrics).to_csv(output_dir / "metrics.csv", index=False)
        pd.DataFrame(predictions).to_csv(output_dir / "predictions_long.csv", index=False)
        pd.DataFrame(uq_rows).to_csv(output_dir / "uq_scores_long.csv", index=False)
        pd.DataFrame(risk_rows).to_csv(output_dir / "risk_coverage.csv", index=False)
        pd.concat(matches, ignore_index=True).to_csv(output_dir / "matching_pairs.csv", index=False)


if __name__ == "__main__":
    main()
