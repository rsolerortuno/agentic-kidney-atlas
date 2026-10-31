"""Command line entry point: `uv run ccrcc <command>`.

run-all rebuilds every deterministic artifact from scratch (data -> omics -> literature tables and validation
-> external enrichment -> graph -> convergence -> assembly -> verification -> HTML -> validation). LLM steps
(relation extraction, agent analyses, narrative) are not re-run: their saved, verified outputs are reused, and
`claude -p "/informe"` (scripts/run_informe.ps1) regenerates them.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time

from .paths import GRAPH_DB, LITERATURE, PROCESSED, RUNS


def _step(name, fn, *a, **kw):
    t0 = time.time()
    print(f"==> {name}", flush=True)
    out = fn(*a, **kw)
    print(f"    done in {time.time() - t0:.0f}s", flush=True)
    return out


def run_all(fresh: bool = True) -> int:
    from . import download, enrich, integrate, kg_build, literature, omics
    from .assemble import assemble
    from .render import render
    from .verify import verify_agents, verify_final

    if fresh:
        shutil.rmtree(PROCESSED, ignore_errors=True)
        PROCESSED.mkdir(parents=True)
        GRAPH_DB.unlink(missing_ok=True)
        for f in ("publications.parquet", "annotations.parquet", "chunks.parquet", "relations.parquet"):
            (LITERATURE / f).unlink(missing_ok=True)
    _step("M1 download (verifies sha256 manifest)", download.download_all)
    _step("M11 download single-cell / spatial h5ad (cached)", download.download_singlecell)
    _step("M2 omics analysis", omics.run)
    _step("M3 literature ingestion (PubTator3, frozen PMID lists)", literature.ingest)
    _step("M3 rebuild extraction batches (deterministic)", literature.write_extraction_batches)
    print("    relations:", _step("M3 validate extracted relations", literature.validate_extractions))
    _step("M5 external enrichment (cached)", enrich.run)
    _step("M4 knowledge graph", kg_build.build)
    from . import kg_cell_structure, singlecell, spatial, structure
    _step("M11 single-cell cell-type expression", singlecell.run)
    _step("M11 spatial deconvolution and neighbourhoods", spatial.run)
    _step("M11 AlphaFold models, pockets and MGD scoring", structure.run)
    print("    claims:", _step("M11 load cell / structure layers into the graph", kg_cell_structure.load)["claims_added"])
    runs = sorted(p.parent for p in RUNS.glob("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]/narrative.json"))
    if not runs:
        print("No saved agent run (runs/<date>/narrative.json). Run scripts/run_informe.ps1 first.")
        return 1
    run = runs[-1]
    _step("convergence", integrate.write, day=run.name)
    print("    agents:", {k: (v.get("kept"), v.get("dropped")) for k, v in
                          _step("verify agents", verify_agents, run)["agents"].items()})
    _step("assemble findings.json", assemble, run)
    rep = _step("verify final", verify_final, run / "findings.json")
    print(f"    statements checked: {rep['statements']} | unsupported claims: {rep['unsupported']}")
    path = _step("M7 render HTML", render, run / "findings.json")
    from importlib import import_module
    sys.path.insert(0, str(PROCESSED.parents[1] / "scripts"))
    code = import_module("validate_report").main(path)
    return 1 if rep["unsupported"] or code else 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="ccrcc")
    sub = p.add_subparsers(dest="cmd", required=True)
    ra = sub.add_parser("run-all", help="rebuild everything deterministic from scratch and render the report")
    ra.add_argument("--keep", action="store_true", help="do not delete processed data / graph first")
    for name in ("download", "omics", "literature", "validate-relations", "enrich", "graph", "integrate", "eval",
                 "singlecell", "spatial", "structure", "cell-graph"):
        sub.add_parser(name)
    ev = sub.choices["eval"]
    ev.add_argument("--qa", action="store_true", help="also run the known-answer questions with claude -p")
    ev.add_argument("--answers", help="grade a saved answers JSON instead of calling claude -p")
    ev.add_argument("--answers-model", default="saved answers")
    args = p.parse_args(argv)

    if args.cmd == "run-all":
        return run_all(fresh=not args.keep)
    if args.cmd == "download":
        from .download import download_all
        print(json.dumps(download_all(), indent=1))
    elif args.cmd == "omics":
        from .omics import run
        print(json.dumps(run(), indent=1, default=str))
    elif args.cmd == "literature":
        from .literature import ingest, write_extraction_batches
        print(json.dumps(ingest(), indent=1)); print(write_extraction_batches())
    elif args.cmd == "validate-relations":
        from .literature import validate_extractions
        print(json.dumps(validate_extractions(), indent=1))
    elif args.cmd == "enrich":
        from .enrich import run
        print(json.dumps(run(), indent=1))
    elif args.cmd == "graph":
        from .kg_build import build
        print(json.dumps(build(), indent=1))
    elif args.cmd == "singlecell":
        from .singlecell import run
        print(json.dumps(run(), indent=1))
    elif args.cmd == "spatial":
        from .spatial import run
        s = run()
        print(json.dumps({k: v for k, v in s.items() if k != "marker_genes"}, indent=1))
    elif args.cmd == "structure":
        from .structure import run
        print(json.dumps(run(), indent=1))
    elif args.cmd == "cell-graph":
        from .kg_cell_structure import load as load_cs
        print(json.dumps(load_cs(), indent=1))
    elif args.cmd == "integrate":
        from .integrate import write
        print(json.dumps(write(), indent=1))
    elif args.cmd == "eval":
        from .evaluate import grade_answers, recall, run_qa, write_results
        rec = recall()
        qa = run_qa() if args.qa else (grade_answers(args.answers, args.answers_model) if args.answers else None)
        print(write_results(rec, qa))
        print(f"recall@{rec['top_n']}: {rec['recall']:.0%}" + (f" | QA {qa['correct']}/{qa['total']}" if qa else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
