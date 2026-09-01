"""M6 NS (case/45) の陰解法 cfl_pseudo 上限探索スイープ。

収束済みの warm 場 (run_0007 [TP] / run_0011 [CPG]) から同一メッシュ restart し、
cfl (= cfl_pseudo) と関連ノブを振って 2000 step 走らせ、発散の有無・初 NaN step・
残差挙動を記録する。probe は run_0015_cflsweep/<tag>/ に置く。

usage: design/.venv-opt/bin/python sweep_cfl_implicit.py <probes.json>
  probes.json: [{"tag": "tp_cfl4", "base": "run_0007_ns_trim_v3", "cfl": 4.0,
                 "patch": {"lowMachPrecond": "2"}}, ...]
  patch は solverConfig の deltaT ブロック内 key: value の置換/追記 (文字列)。
出力: run_0015_cflsweep/summary.json (+ 各 probe の residual_history / res)
"""
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate.runner import FORGE_TOOLS, _ENV, run_forge  # noqa: E402

CASE = Path(__file__).resolve().parent
SWEEP = CASE / "run_0015_cflsweep"
NSTEP = 2000
OUTINT = 500


def make_probe(tag: str, base: str, cfl: float, patch: dict | None,
               ic_from: str | None = None, nstep: int = NSTEP) -> Path:
    rd = SWEEP / tag
    if rd.exists():
        shutil.rmtree(rd)
    rd.mkdir(parents=True)
    b = CASE / base
    for f in ("solverConfig.yaml", "bcondConfig.yaml", "probe.yaml",
              "species_db.yaml", "nozzle.h5"):
        if (b / f).exists():
            shutil.copy(b / f, rd / f)
    # 収束場を restart (ic_from 指定時は cross-mesh interp)
    src = CASE / ic_from if ic_from else b
    res = sorted(src.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                    str(res), str(rd / "nozzle.h5")],
                   env=_ENV, check=True, capture_output=True, text=True)
    cfg = (rd / "solverConfig.yaml").read_text()
    cfg = re.sub(r"nStepOuter: \d+", f"nStepOuter: {nstep}", cfg)
    cfg = re.sub(r"outStepInterval: \d+", f"outStepInterval: {OUTINT}", cfg)
    cfg = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+",
                 f"cfl: {cfl}, cfl_pseudo: {cfl}", cfg)
    for k, v in (patch or {}).items():
        if re.search(rf"{k}: [^,}}\n]+", cfg):
            cfg = re.sub(rf"{k}: [^,}}\n]+", f"{k}: {v}", cfg)
        else:  # deltaT ブロックに追記
            cfg = cfg.replace("blockDPLUR: 1", f"blockDPLUR: 1, {k}: {v}", 1)
    (rd / "solverConfig.yaml").write_text(cfg)
    return rd


def analyze(rd: Path) -> dict:
    out = {"finished": False, "first_nan_step": None, "max_growth": None,
           "final": {}}
    f = rd / "residual_history.csv"
    if not f.exists():
        out["error"] = "no residual_history"
        return out
    d = np.genfromtxt(f, delimiter=",", names=True)
    if d.size == 0:
        out["error"] = "empty history"
        return out
    step = d["step"] if "step" in d.dtype.names else np.arange(len(d))
    for n in d.dtype.names:
        if not n.startswith("rms"):
            continue
        v = np.atleast_1d(d[n])
        bad = ~np.isfinite(v)
        if bad.any():
            s = int(np.atleast_1d(step)[np.argmax(bad)])
            if out["first_nan_step"] is None or s < out["first_nan_step"]:
                out["first_nan_step"] = s
                out["first_nan_field"] = n
        else:
            out["final"][n] = float(v[-1])
            if v[0] > 0:
                g = float(np.max(v) / v[0])
                if out["max_growth"] is None or g > out["max_growth"]:
                    out["max_growth"] = g
                    out["max_growth_field"] = n
    out["finished"] = out["first_nan_step"] is None
    out["last_step"] = int(np.atleast_1d(step)[-1])
    return out


def main():
    probes = json.loads(Path(sys.argv[1]).read_text())
    SWEEP.mkdir(exist_ok=True)
    summ_path = SWEEP / "summary.json"
    summary = json.loads(summ_path.read_text()) if summ_path.exists() else {}
    for p in probes:
        tag = p["tag"]
        print(f"=== {tag}: base={p['base']} cfl={p['cfl']} patch={p.get('patch')}",
              flush=True)
        rd = make_probe(tag, p["base"], p["cfl"], p.get("patch"),
                        p.get("ic_from"), int(p.get("nstep", NSTEP)))
        rc = run_forge(rd)
        a = analyze(rd)
        a.update({"rc": rc, "cfl": p["cfl"], "base": p["base"],
                  "patch": p.get("patch")})
        summary[tag] = a
        verdict = "OK(完走)" if a["finished"] else \
            f"DIVERGED@{a.get('first_nan_step')} ({a.get('first_nan_field')})" \
            if a.get("first_nan_step") is not None else f"STOP@{a.get('last_step')}"
        print(f"    -> {verdict}  rc={rc}  max_growth={a.get('max_growth'):.3g} "
              f"({a.get('max_growth_field')})", flush=True)
        # 大きい res は掃除 (最後の 1 つと発散直前だけ残す)
        res = sorted(rd.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        for f in res[:-2]:
            f.unlink()
        summ_path.write_text(json.dumps(summary, indent=1))
    print("saved:", summ_path)


if __name__ == "__main__":
    main()
