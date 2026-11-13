"""Run the MoRank manuscript experiments from the repository root."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
SIM = ROOT / "simulation"
REAL = ROOT / "real_data"
ABLATIONS = ("DropOriginalWSKR", "DropMomentBlock", "DropRankBlock", "RankOnly")
REAL_STAGES = ("discovery", "bulk", "singlecell", "string", "context")


def verify_models():
    for line in (SIM / "model_sha256.txt").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split(maxsplit=1)
        digest = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if digest != expected:
            raise ValueError(f"Frozen model checksum mismatch: {relative}")


class Runner:
    def __init__(self, args):
        self.output = args.output.resolve()
        if self.output == ROOT or any(
            self.output == p or p in self.output.parents for p in (SIM, REAL)
        ):
            raise ValueError("Use an output directory outside the experiment source folders.")
        self.workers = args.workers
        self.env = os.environ.copy()
        self.env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        if os.name == "nt":
            self.env.update(LC_ALL=".UTF-8", LC_CTYPE=".UTF-8")
        self.rscript = args.rscript
        self.primary = self.output / "primary"
        self.common = self.output / "common"
        self.replication = self.output / "replication"
        self.robustness = self.output / "robustness"
        self.null = self.output / "null"
        self.env.update(
            MORANK_REPO_DIR=str(ROOT),
            MORANK_PRIMARY_DIR=str(self.primary),
            MORANK_COMMON_DIR=str(self.common),
            MORANK_AUDIT_REPLICATION=str(self.replication),
            MORANK_DEVELOPMENT=str(self.output / "training" / "reference_features"),
            MORANK_MIXED_CACHE=str(self.output / "training" / "mixed_features"),
        )
        for path in (
            self.output / "results" / "final_model", self.output / "qa",
            self.output / "logs", self.common / "results", self.robustness,
            self.replication / "models", self.replication / "results", self.null,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def command(self, executable, directory, script, *arguments, **extra_env):
        details = [str(arg) for arg in arguments
                   if isinstance(arg, int) or isinstance(arg, str) and arg.replace("_", "").isalnum()]
        label = "-".join([Path(script).stem, *details])
        log = self.output / "logs" / f"{label}.log"
        cmd = [*executable, script, *map(str, arguments)]
        with log.open("w", encoding="utf-8") as stream:
            result = subprocess.run(
                cmd, cwd=directory, env={**self.env, **{k: str(v) for k, v in extra_env.items()}},
                stdout=stream, stderr=subprocess.STDOUT, check=False,
            )
        if result.returncode:
            raise RuntimeError(f"{script} failed (exit {result.returncode}). See {log}")
        print(f"Completed {label}", flush=True)

    def r(self, directory, script, *arguments, **env):
        self.command([self.rscript, "--vanilla"], directory, script, *arguments, **env)

    def py(self, directory, script, *arguments):
        self.command([sys.executable, "-B"], directory, script, *arguments)

    def batch(self, count, function):
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            list(pool.map(function, range(1, count + 1)))

    def smoke(self):
        self.r(SIM / "reference", "test_redesign.R")
        self.r(SIM / "reference", "test_uniform_permutation.R")
        self.r(SIM / "common", "check_formulas.R", self.output)
        self.r(SIM / "robustness", "smoke_frozen.R")
        self.r(SIM / "reference", "run_final_benchmark.R", 1, 24, self.output / "smoke", 1)

    def simulation(self, stage):
        ref, rob, common = SIM / "reference", SIM / "robustness", SIM / "common"
        if stage in ("primary", "all"):
            self.batch(12, lambda i: self.r(ref, "run_final_benchmark.R", i, 12, self.primary))
            self.r(ref, "summarize_final.R", self.primary)
        if stage in ("ablations", "all"):
            self.r(SIM / "ablations", "evaluate_ablation.R",
                   MORANK_AUDIT_RESULTS=self.primary / "results")
            self.r(SIM / "ablations", "summarize.R", self.primary)
        if stage in ("stress", "all"):
            self.r(ref, "run_stress.R", self.primary)
            self.batch(12, lambda i: self.r(rob, "evaluate_final_figure.R", "replay", i, 12, self.null))
            self.r(rob, "summarize_stress_and_null.R", self.output, "stress")
        if stage in ("robustness", "all"):
            self.batch(12, lambda i: self.r(rob, "evaluate_robustness.R", "test2", i, 12,
                                            MORANK_AUDIT_RESULTS=self.robustness))
            self.r(rob, "summarize_test2_and_common.R", self.output, "test2")
        if stage in ("replication", "all"):
            shutil.copyfile(SIM / "replication" / "selected_gate.csv",
                            self.replication / "models" / "selected_gate.csv")
            self.batch(20, lambda i: self.r(rob, "evaluate_replication.R", "evaluate_new", i, 20))
            self.r(self.replication, str(SIM / "replication" / "summarize.R"))
        if stage == "routing":
            for split in ("profile_fit", "profile_tune"):
                self.batch(6, lambda i: self.r(rob, "evaluate_replication.R", split, i, 6))
            self.r(rob, "evaluate_replication.R", "select_gate")
        if stage in ("null", "all"):
            self.batch(26, lambda i: self.r(rob, "evaluate_final_figure.R", "calibration", i, 26, self.null))
            self.r(rob, "summarize_stress_and_null.R", self.output, "null")
        if stage in ("common", "all"):
            self.batch(6, lambda i: self.r(common, "prepare_common_genes.R", i))
            self.batch(24, lambda i: self.r(common, "evaluate_waddr.R", i))
            self.batch(24, lambda i: self.py(common, "evaluate_mmdagg.py", i))
            self.batch(6, lambda i: self.r(rob, "evaluate_robustness.R", "neighbor_final", i, 6,
                                           MORANK_AUDIT_RESULTS=self.common / "results"))
            self.r(rob, "summarize_test2_and_common.R", self.output, "neighbor_final")

    def train(self):
        ref, rob = SIM / "reference", SIM / "robustness"
        output = self.output / "training"
        branch = output / "branches"
        gate = output / "gate"
        for folder in (branch, gate):
            for child in ("models", "results"):
                (folder / child).mkdir(parents=True, exist_ok=True)
        for split in ("fit", "validation"):
            self.r(ref, "generate_development.R", split)
        for family in ("ClassicalStack", "WSKR_A"):
            self.r(ref, "fit_reference.R", family, MORANK_TRAIN_OUTPUT=output / "reference")
        for family in ABLATIONS:
            self.r(SIM / "ablations", "fit_ablations.R", family, MORANK_TRAIN_OUTPUT=output / "ablations")
        for split in ("fit", "tune"):
            self.batch(8, lambda i: self.r(rob, "generate.R", split, i, 8))
        for family in ("Classical", "MoRank37", "MoRank48"):
            self.r(rob, "audit_mixed_fit.R", "train", family, MORANK_FIT_OUTPUT=branch)
        self.r(rob, "audit_mixed_fit.R", "select", MORANK_FIT_OUTPUT=branch)
        for family in ("Classical", "MoRank37"):
            self.r(rob, "audit_specialist_fit.R", "specialist", family, MORANK_FIT_OUTPUT=branch)
        self.r(rob, "audit_initial_gate.R", "train_gate", MORANK_AUDIT_GATE_OUTPUT=gate)

    def real(self, args):
        data_root = args.input_dir.resolve()
        if not data_root.is_dir():
            raise ValueError(f"Input directory does not exist: {data_root}")
        self.env["MORANK_INPUT_DIR"] = str(data_root)
        tcga = args.tcga_rdata.resolve() if args.tcga_rdata else data_root / "TCGA-LUAD_list.RData"
        self.env["MORANK_TCGA_RDATA"] = str(tcga)
        stages = REAL_STAGES if args.stage == "all" else (args.stage,)
        for stage in stages:
            if stage == "string":
                self.py(REAL, "final_string.py", self.output)
            elif stage == "ko":
                if not args.alignment_r or not args.alignment_r.is_file():
                    raise ValueError("The ko stage requires --alignment-r with the reference alignment function.")
                self.env["MORANK_ALIGNMENT_R"] = str(args.alignment_r.resolve())
                self.batch(8, lambda i: self.r(REAL, "audit_ko_alignment.R", self.output, i, 8))
            elif stage == "singlecell-input":
                self.r(REAL, "audit_singlecell_input.R", self.output)
            elif stage in ("tcga_deseq2", "tcga_edger", "tcga_aux", "bulk-comparators"):
                mode = "bulk" if stage == "bulk-comparators" else stage
                self.r(REAL, "audit_count_features.R", self.output, mode)
            else:
                self.r(REAL, f"final_{stage}.R", self.output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", "Rscript"))
    parser.add_argument("--output", type=Path, default=ROOT / "runs")
    parser.add_argument("--workers", type=int, default=2, help="Maximum concurrent workers")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("smoke")
    commands.add_parser("train")
    simulation = commands.add_parser("simulation")
    simulation.add_argument("stage", choices=(
        "all", "primary", "ablations", "stress", "robustness", "replication", "null", "common", "routing"))
    real = commands.add_parser("real")
    real.add_argument("stage", choices=(
        "all", *REAL_STAGES, "ko", "singlecell-input", "tcga_deseq2",
        "tcga_edger", "tcga_aux", "bulk-comparators"))
    real.add_argument("--input-dir", type=Path, required=True)
    real.add_argument("--tcga-rdata", type=Path)
    real.add_argument("--alignment-r", type=Path)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    try:
        verify_models()
        runner = Runner(args)
        if args.command == "simulation":
            runner.simulation(args.stage)
        elif args.command == "real":
            runner.real(args)
        else:
            getattr(runner, args.command)()
    except (OSError, RuntimeError, ValueError) as exc:
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
