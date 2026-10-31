"""Command line entry point: `uv run ccrcc <command>`.

run-all rebuilds every deterministic artifact from scratch (data -> omics -> clinical/HPA -> M13 extra layers ->
literature tables and validation -> external enrichment -> graph -> cell/structure/discovery layers -> convergence v3,
profiles and sensitivity -> assembly -> semantic verification -> HTML -> validation). LLM steps
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
    from . import extra_download, kg_layers, layers, safety_clinical
    _step("M1 download (verifies sha256 manifest)", download.download_all)
    _step("M11 download single-cell / spatial h5ad (cached)", download.download_singlecell)
    _step("M13 download extra layers (mutations, DepMap, TCGA, phosphosites, purity, OT bulk; cached)",
          extra_download.download_extra)
    _step("M2 omics analysis", omics.run)
    # HPA normal tissue + CPTAC grade/stage association: needed by M12 discovery and the v3 biomarker profile.
    # (Before M13 run-all never ran this step, so a clean run failed at M12 unless the files were already there.)
    _step("M12 HPA vital-organ expression + clinical association", safety_clinical.run)
    print("    layers:", {k: v for k, v in _step("M13 extra layers (mutations, cis dosage, DepMap, TCGA, kinases, purity)",
                                               layers.run).items() if k in ("mutations", "tcga", "purity")})
    _step("M3 literature ingestion (PubTator3, frozen PMID lists)", literature.ingest)
    _step("M3 rebuild extraction batches (deterministic)", literature.write_extraction_batches)
    print("    relations:", _step("M3 validate extracted relations", literature.validate_extractions))
    _step("M5 external enrichment (cached)", enrich.run)
    _step("M4 knowledge graph", kg_build.build)
    print("    M13 claims:", _step("M13 load extra layers into the graph", kg_layers.load_layers)["claims_added"])
    from . import discovery, kg_cell_structure, singlecell, spatial, structure
    _step("M11 single-cell cell-type expression", singlecell.run)
    print("    discovery:", {k: v for k, v in _step("M12 target discovery (Europe PMC, Open Targets, UniProt)",
                                                   discovery.run).items() if k in ("pool", "novel", "by_status")})
    _step("M12 single-cell summary for report + discovery shortlist", singlecell.write_summary)
    _step("M11 spatial deconvolution and neighbourhoods", spatial.run)
    _step("M12 spatial section maps for the report", spatial.section_maps)
    _step("M11 AlphaFold models, pockets and MGD scoring", structure.run)
    print("    claims:", _step("M11 load cell / structure layers into the graph", kg_cell_structure.load)["claims_added"])
    print("    discovery claims:", _step("M12 load discovery profiles into the graph", discovery.load_graph))
    runs = sorted(p.parent for p in RUNS.glob("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]/narrative.json"))
    if not runs:
        print("No saved agent run (runs/<date>/narrative.json). Run scripts/run_informe.ps1 first.")
        return 1
    run = runs[-1]
    _step("convergence v3 + target / biomarker profiles + weight sensitivity", integrate.write, day=run.name)
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
                 "singlecell", "spatial", "structure", "cell-graph", "discovery", "layers", "layers-graph",
                 "download-extra", "safety-clinical", "baselines"):
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
    elif args.cmd == "discovery":
        from .discovery import run
        print(json.dumps(run(), indent=1))
    elif args.cmd == "cell-graph":
        from .kg_cell_structure import load as load_cs
        print(json.dumps(load_cs(), indent=1))
    elif args.cmd == "download-extra":
        from .extra_download import download_extra
        print(json.dumps(download_extra(), indent=1))
    elif args.cmd == "safety-clinical":
        from .safety_clinical import run
        print(json.dumps(run(), indent=1))
    elif args.cmd == "layers":
        from .layers import run
        print(json.dumps(run(), indent=1, default=float))
    elif args.cmd == "layers-graph":
        from .kg_layers import load_layers
        print(json.dumps(load_layers(), indent=1))
    elif args.cmd == "baselines":
        from .baselines import run
        print(json.dumps(run(), indent=1, default=str))
    elif args.cmd == "integrate":
        from .integrate import write
        print(json.dumps(write(), indent=1))
    elif args.cmd == "eval":
        from .evaluate import EVALS, grade_answers, recall, run_qa, write_results
        from .integrate import compute
        from .outcomes import run as outcome_eval
        rec = recall()
        qa = run_qa() if args.qa else (grade_answers(args.answers, args.answers_model) if args.answers else None)
        if qa is None and (EVALS / "results.json").exists():   # keep the last QA run unless asked to re-run it
            qa = json.loads((EVALS / "results.json").read_text(encoding="utf-8")).get("qa")
        outc = outcome_eval(compute())
        sens_p = PROCESSED / "sensitivity_summary.json"
        sens = json.loads(sens_p.read_text()) if sens_p.exists() else None
        base_p = EVALS / "baselines" / "baselines.json"
        base = json.loads(base_p.read_text(encoding="utf-8")) if base_p.exists() else None
        print(write_results(rec, qa, outc, base, sens))
        print(f"recall@{rec['top_n']}: {rec['recall']:.0%}" + (f" | QA {qa['correct']}/{qa['total']}" if qa else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
