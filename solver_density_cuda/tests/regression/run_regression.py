#!/usr/bin/env python3
"""forge 回帰テストハーネス.

短時間で回せる検証ケースを 1 本走らせ、生成された ``residual_history.csv`` を
コミット済みの baseline と相対許容差で比較して回帰の有無を判定する。

GPU 必須。実行は AGENTS.md の共通ルールに従い、既存 run を使い回さず複製した
新しい run ディレクトリで行う (本ハーネスがケース内に一時 run を作って実行する)。

モード:
  (既定)            baseline と比較して PASS/FAIL を返す
  --case all        cases/*.json を名前順にすべて実行する
  --smoke           短ステップで実行し、残差の有限性のみ確認 (baseline 不要)
  --update-baseline 実行結果を baseline として保存し直す (意図的変更後に使用)
  --compare-only F  実行せず、既存の residual_history.csv (F) を baseline と比較

runner:
  --runner auto|native|docker  (既定 auto: native バイナリがあれば native)

終了コード: 0=PASS, 1=FAIL(回帰検出), 2=セットアップ/実行エラー
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# .../solver_density_cuda/tests/regression/run_regression.py
THIS_DIR = Path(__file__).resolve().parent
SOLVER_DIR = THIS_DIR.parents[1]          # solver_density_cuda
REPO_ROOT = THIS_DIR.parents[2]           # リポジトリルート
CASES_DIR = THIS_DIR / "cases"
BASELINES_DIR = THIS_DIR / "baselines"

# native forge バイナリ候補 (優先順)
NATIVE_FORGE_CANDIDATES = [
    SOLVER_DIR / ".build-native" / "relwithdebinfo" / "forge",
    SOLVER_DIR / "build" / "forge",
]
DOCKER_FORGE = "/workspace/solver_density_cuda/build/forge"
DOCKER_IMAGE_DEFAULT = "forge-solver:cuda-dev"
WSL_CUDA_LIB_DIR = "/usr/lib/wsl/lib"


# --------------------------------------------------------------------------- #
# ケース設定
# --------------------------------------------------------------------------- #
def load_case(name: str) -> dict:
    path = CASES_DIR / f"{name}.json"
    if not path.exists():
        die(f"ケース設定が見つからない: {path}")
    with path.open(encoding="utf-8") as f:
        cfg = json.load(f)
    cfg["_path"] = str(path)
    return cfg


def list_case_names() -> list[str]:
    return sorted(p.stem for p in CASES_DIR.glob("*.json"))


# --------------------------------------------------------------------------- #
# 残差 CSV の読み込み・比較
# --------------------------------------------------------------------------- #
def load_outer_rows(csv_path: Path) -> dict:
    """outer_begin 行 (inner == -1) を step -> {col: value} の辞書で返す."""
    data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=float,
                         encoding="utf-8")
    if data.ndim == 0:
        data = data.reshape(1)
    mask = data["inner"] == -1
    rows = data[mask]
    out = {}
    for r in rows:
        out[int(r["step"])] = r
    return out


def check_finite(rows: dict, columns: list) -> list:
    """全 step・全比較列が有限か確認し、違反を文字列リストで返す."""
    problems = []
    for step in sorted(rows):
        for col in columns:
            v = rows[step][col]
            if not np.isfinite(v):
                problems.append(f"step {step}: {col} = {v} (非有限)")
    return problems


def compare(baseline: dict, candidate: dict, columns: list,
            rtol: float, atol: float, sample_steps: list) -> list:
    """baseline と candidate を比較し、許容差超過を文字列リストで返す."""
    failures = []
    final_step = max(baseline)
    steps_to_check = sorted(set([s for s in sample_steps if s in baseline] + [final_step]))
    for step in steps_to_check:
        if step not in candidate:
            failures.append(f"step {step}: candidate に該当行がない")
            continue
        b = baseline[step]
        c = candidate[step]
        for col in columns:
            bv = float(b[col])
            cv = float(c[col])
            if not np.isfinite(cv):
                failures.append(f"step {step} {col}: candidate が非有限 ({cv})")
                continue
            tol = atol + rtol * abs(bv)
            diff = abs(cv - bv)
            if diff > tol:
                rel = diff / abs(bv) if bv != 0 else float("inf")
                failures.append(
                    f"step {step} {col}: baseline={bv:.6e} candidate={cv:.6e} "
                    f"diff={diff:.3e} (許容={tol:.3e}, rel={rel:.2%})")
    return failures


# --------------------------------------------------------------------------- #
# run ディレクトリの準備と forge 実行
# --------------------------------------------------------------------------- #
CONFIG_FILES = ["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml"]


def prepare_run_dir(case_dir: Path, template_run: str, mesh_file: str,
                    tag: str, nstep: int) -> Path:
    """雛型 run を複製した新しい run ディレクトリを作り、nStep を書き換える."""
    template = case_dir / template_run
    if not template.is_dir():
        die(f"雛型 run が見つからない: {template}")
    mesh_path = template / mesh_file
    if not mesh_path.exists():
        die(f"メッシュ {mesh_file} が雛型 run に無い。先に run_case.sh 等でメッシュを生成すること: {mesh_path}")

    run_dir = case_dir / f"run_regression_{tag}"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)

    # 必要最小限のファイルのみ複製 (res_*.h5 やログは持ち込まない)
    for name in CONFIG_FILES:
        src = template / name
        if src.exists():
            shutil.copy2(src, run_dir / name)
    shutil.copy2(mesh_path, run_dir / mesh_file)

    # valueFileName が mesh_file と別名のケース (TP など) は元ファイルも複製する。
    solver_config_src = template / "solverConfig.yaml"
    if solver_config_src.exists():
        value_file = _extract_value_file_name(solver_config_src)
        if value_file and value_file != mesh_file:
            value_src = template / value_file
            if value_src.exists():
                shutil.copy2(value_src, run_dir / value_file)

    prepare_solver_config(run_dir / "solverConfig.yaml", nstep)
    return run_dir


def _extract_value_file_name(solver_config: Path) -> str | None:
    """solverConfig.yaml から valueFileName を雑に抜き出す."""
    text = solver_config.read_text(encoding="utf-8")
    m = re.search(r'^\s*valueFileName\s*:\s*"?([^"\n#]+)"?\s*$', text, re.MULTILINE)
    if not m:
        return None
    return m.group(1).strip()


def prepare_solver_config(solver_config: Path, nstep: int) -> None:
    """旧キー名を現行ソルバの名称へ正規化し、外側ステップ数を nstep に設定する.

    現行 forge は ``nStep`` → ``nStepOuter``、``nInnerLoop`` / ``dualTime_InnerLoop``
    → ``nStepInner`` を要求する (input/solverConfig.cpp 参照)。既存ケース設定が
    旧キーのままでも回帰テストを回せるよう、コメントを保持したまま置換する。
    """
    text = solver_config.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^(\s*)nInnerLoop(\s*:)", r"\1nStepInner\2", text)
    text = re.sub(r"(?m)^(\s*)dualTime_InnerLoop(\s*:)", r"\1nStepInner\2", text)
    text = re.sub(r"(?m)^(\s*)nStep(\s*:)", r"\1nStepOuter\2", text)
    text, n = re.subn(r"(?m)^(\s*nStepOuter\s*:\s*)\d+",
                      lambda m: f"{m.group(1)}{nstep}", text, count=1)
    if n == 0:
        die(f"solverConfig.yaml に nStep/nStepOuter 行が見つからない: {solver_config}")
    solver_config.write_text(text, encoding="utf-8")


def resolve_native_forge() -> Path | None:
    for cand in NATIVE_FORGE_CANDIDATES:
        if cand.exists():
            return cand
    return None


def run_forge(run_dir: Path, runner: str, case_rel: str, image: str) -> None:
    """run_dir 内で forge を実行する."""
    if runner == "native":
        binary = resolve_native_forge()
        if binary is None:
            die("native forge バイナリが見つからない。tools/build_native_wsl.sh でビルドすること")
        cmd = [str(binary)]
        print(f"[run] native: {binary} (cwd={run_dir})")
        rc = subprocess.run(cmd, cwd=run_dir, env=native_run_env()).returncode
    elif runner == "docker":
        workdir = f"/workspace/{case_rel}/{run_dir.name}"
        cmd = [
            "docker", "run", "--rm", "--gpus", "all",
            "-v", f"{REPO_ROOT}:/workspace",
            "-w", workdir,
            image, "bash", "-lc", DOCKER_FORGE,
        ]
        print(f"[run] docker: {' '.join(cmd)}")
        rc = subprocess.run(cmd).returncode
    else:
        die(f"未知の runner: {runner}")
    if rc != 0:
        die(f"forge 実行が失敗 (exit={rc})")


def generate_png(run_dir: Path, columns: list) -> None:
    """residual_history.csv から residual_history.png を生成 (失敗は警告のみ).

    AGENTS.md の「残差 PNG を残す」ルールに従う。explicit/implicit で残差 CSV の
    列名が異なるため、外部ツールに依存せず outer_begin 行の比較列を直接描画する。
    """
    csv = run_dir / "residual_history.csv"
    png = run_dir / "residual_history.png"
    if not csv.exists():
        print(f"[warn] {csv} が無いため PNG 生成をスキップ")
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rows = load_outer_rows(csv)
        steps = sorted(rows)
        xs = np.array(steps, dtype=float)
        fig, ax = plt.subplots(figsize=(8, 5))
        for col in columns:
            ys = np.array([rows[s][col] for s in steps], dtype=float)
            ys = np.where(ys > 0, ys, np.nan)  # log 軸用に非正値を除外
            ax.semilogy(xs, ys, label=col)
        ax.set_xlabel("step")
        ax.set_ylabel("RMS residual")
        ax.set_title(f"residual history ({run_dir.name})")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(png, dpi=120)
        plt.close(fig)
        print(f"[ok] PNG 生成: {png}")
    except Exception as exc:  # noqa: BLE001 - PNG 生成失敗はテストを落とさない
        print(f"[warn] PNG 生成に失敗 (テストは継続): {exc}")


# --------------------------------------------------------------------------- #
# ユーティリティ
# --------------------------------------------------------------------------- #
def die(msg: str) -> None:
    print(f"[error] {msg}", file=sys.stderr)
    sys.exit(2)


def decide_runner(requested: str) -> str:
    if requested != "auto":
        return requested
    return "native" if resolve_native_forge() is not None else "docker"


def native_run_env() -> dict[str, str]:
    """WSL native 実行時の CUDA ライブラリ探索順を整える."""
    env = os.environ.copy()
    if os.path.exists(WSL_CUDA_LIB_DIR):
        current = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = (
            f"{WSL_CUDA_LIB_DIR}:{current}" if current else WSL_CUDA_LIB_DIR
        )
    return env


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description="forge 回帰テストハーネス")
    ap.add_argument("--case", default="naca_slau",
                    help="ケース設定名 (cases/<name>.json)。all で全ケース")
    ap.add_argument("--runner", default="auto", choices=["auto", "native", "docker"])
    ap.add_argument("--image", default=DOCKER_IMAGE_DEFAULT, help="docker イメージ名")
    ap.add_argument("--smoke", action="store_true", help="短ステップで有限性のみ確認")
    ap.add_argument("--update-baseline", action="store_true", help="実行結果を baseline に保存")
    ap.add_argument("--compare-only", metavar="CSV",
                    help="実行せず、指定 CSV を baseline と比較")
    ap.add_argument("--keep", action="store_true", help="一時 run ディレクトリを残す")
    ap.add_argument("--tag", default=None, help="一時 run ディレクトリ名のサフィックス")
    args = ap.parse_args()

    if args.case == "all":
        if args.compare_only:
            die("--compare-only は --case all と併用できない")
        failures = []
        for name in list_case_names():
            print(f"\n===== case: {name} =====")
            try:
                rc = run_one_case(name, args)
            except SystemExit as exc:
                rc = int(exc.code) if isinstance(exc.code, int) else 2
            if rc != 0:
                failures.append(name)
        if failures:
            print(f"\n[FAIL] all: 失敗ケース {', '.join(failures)}")
            return 1
        print("\n[PASS] all: 全ケース問題なし")
        return 0

    return run_one_case(args.case, args)


def run_one_case(case_name: str, args: argparse.Namespace) -> int:
    """1 ケース分を実行する。--case all からも再利用する。"""
    cfg = load_case(case_name)
    columns = cfg["compare_columns"]
    baseline_csv = BASELINES_DIR / cfg["name"] / "residual_history.csv"

    # --- compare-only: GPU 不要、比較ロジックのみ ---
    if args.compare_only:
        if not baseline_csv.exists():
            die(f"baseline が無い: {baseline_csv}")
        baseline = load_outer_rows(baseline_csv)
        candidate = load_outer_rows(Path(args.compare_only))
        failures = compare(baseline, candidate, columns,
                           cfg["rtol"], cfg["atol"], cfg.get("sample_steps", []))
        return report(failures)

    # --- 実行系 (GPU 必須) ---
    runner = decide_runner(args.runner)
    case_dir = REPO_ROOT / cfg["case_dir"]
    tag_base = args.tag or ("smoke" if args.smoke else "regression")
    tag = tag_base if case_name == args.case else f"{tag_base}_{case_name}"
    nstep = cfg["nstep_smoke"] if args.smoke else cfg["nstep_regression"]

    run_dir = prepare_run_dir(case_dir, cfg["template_run"], cfg["mesh_file"], tag, nstep)
    try:
        run_forge(run_dir, runner, cfg["case_dir"], args.image)
        generate_png(run_dir, columns)
        produced = run_dir / "residual_history.csv"
        if not produced.exists():
            die(f"residual_history.csv が生成されなかった: {produced}")

        if args.update_baseline:
            baseline_csv.parent.mkdir(parents=True, exist_ok=True)
            _write_baseline(produced, baseline_csv, nstep)
            print(f"[ok] baseline を更新: {baseline_csv}")
            return 0

        candidate = load_outer_rows(produced)

        if args.smoke:
            problems = check_finite(candidate, columns)
            final = max(candidate)
            print(f"[smoke] final step={final} "
                  + " ".join(f"{c}={candidate[final][c]:.3e}" for c in columns))
            return report(problems, label="smoke")

        if not baseline_csv.exists():
            die(f"baseline が無い: {baseline_csv} (先に --update-baseline で生成すること)")
        baseline = load_outer_rows(baseline_csv)
        failures = compare(baseline, candidate, columns,
                           cfg["rtol"], cfg["atol"], cfg.get("sample_steps", []))
        return report(failures)
    finally:
        if not args.keep and run_dir.exists():
            shutil.rmtree(run_dir)
            print(f"[cleanup] {run_dir} を削除 (残すには --keep)")


def _write_baseline(produced: Path, baseline_csv: Path, nstep: int) -> None:
    """生成 CSV から outer_begin 行 (step<=nstep) を抜き出して baseline に保存."""
    lines = produced.read_text(encoding="utf-8").splitlines()
    header = lines[0]
    kept = [header]
    for line in lines[1:]:
        parts = line.split(",")
        if len(parts) < 3:
            continue
        if parts[2] == "outer_begin" and int(float(parts[0])) <= nstep:
            kept.append(line)
    baseline_csv.write_text("\n".join(kept) + "\n", encoding="utf-8")


def report(problems: list, label: str = "regression") -> int:
    if problems:
        print(f"\n[FAIL] {label}: {len(problems)} 件の問題")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"\n[PASS] {label}: 問題なし")
    return 0


if __name__ == "__main__":
    sys.exit(main())
