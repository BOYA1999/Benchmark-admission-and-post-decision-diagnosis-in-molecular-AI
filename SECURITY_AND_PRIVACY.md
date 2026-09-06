# Repository privacy and data handling

This repository is assembled from a fixed public-release allowlist. Its verification checks cover local absolute paths, personal contact strings, common credential formats, unsafe paths, generated caches, checkpoints, model-weight formats, structure-bearing table columns, and document artifacts. Compressed CSV contents are scanned as well. Third-party names and URLs remain only where needed for scientific attribution or licence compliance.

The distributed scientific records concern molecular benchmark data. Structure-bearing upstream tables, individual activities and predictions, local execution state, and weight-dependent caches are kept outside the repository. Public aggregate tables and integer-index matching records retain provenance metadata. Source-level reconstructions belong in a separate working copy, not the upload folder. No directory contains more than 100 directly contained files.

Run `python -B verify_public_package.py` after any change and before creating a release archive.
