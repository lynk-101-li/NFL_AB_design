# Real model runs

This directory is the durable boundary for reviewed RFantibody, IgGM, Germinal,
and downstream structure-model handoffs and results. The deterministic proxy
workflow writes reproducible tables under `outputs/` and must never remove files
from this directory.

Recommended layout:

```text
real_runs/
  handoffs/<source-run-id>/<conformation-id>/<profile>/
  results/<source-run-id>/<engine-or-tool>/<job-id>/
  final/<source-run-id>/
```

Every real run should retain its target-structure manifest, normalized request
hashes, exact argv, upstream revisions, checkpoint hashes, logs, status, and
selected results. Large server-side runs remain under the matching AutoDL topic
directory and may be copied here only when intentionally archived.
