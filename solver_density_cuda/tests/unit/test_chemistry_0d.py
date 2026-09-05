#!/usr/bin/env python3
"""0-D 化学 (有限速度反応) のホスト回帰テスト: forge の反応表/Jacobian/反応器 (tools/test_chemistry.cpp) を
Cantera 参照 (tools/chem_reference_cantera.py) と着火遅れで照合する。GPU 不要。

  python3 tests/unit/test_chemistry_0d.py            # ビルド → 各機構で照合 → PASS/FAIL (exit 0/1)
  python3 tests/unit/test_chemistry_0d.py --update   # Cantera 参照 (goldens/chem0d.json) を再生成 (.venv-chem 必須)

判定: 量論 H₂-air 1200 K 1 atm 定積反応器の着火遅れ (max dT/dt) が Cantera と rtol 以内、Jacobian の有限差分照合 PASS、
Σω≈0。BDF1 の刻みは TCHEM_* で締める (既定刻みは着火遅れを 5–15 % 短く出す: plan chemistry §9 2026-09-05 (7)(e))。
参照値 (Cantera の着火遅れ・終端温度) は committed golden tests/unit/goldens/chem0d.json (機構名 → {tau_ign_s, T_end})。
機構を追加したら MECHS に足して --update する (.venv-chem の Cantera で再生成)。
"""
import argparse, os, re, subprocess, sys, pathlib, shutil

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]                       # solver_density_cuda/
MECH_DIR = ROOT / "tools/mechanisms"
GOLD = HERE / "goldens"
VENV_PY = ROOT.parents[0] / ".venv-chem/bin/python"      # forge-chem/.venv-chem (無ければ ../forge/.venv-chem)
if not VENV_PY.exists():
    VENV_PY = pathlib.Path("/home/sano/work/forge/.venv-chem/bin/python")
MECHS = {  # name: (mechanism yaml, rtol on tau_ign)
    "jachimowski9": (MECH_DIR / "h2o2_jachimowski1988_9sp20r.yaml", 0.02),
    "li2004":       (MECH_DIR / "h2co_li2004_cantera.yaml", 0.02),
}
DB = MECH_DIR / "species_db_h2air_cea.yaml"
T0, P_ATM, T_END = "1200", "1", "5e-4"
ENV = dict(os.environ, TCHEM_RELMAX="0.03", TCHEM_DTMAX="4", TCHEM_DTCAP="5e-8")


def build(bindir: pathlib.Path) -> pathlib.Path:
    exe = bindir / "tchem"
    src = ROOT / "tools/test_chemistry.cpp"
    if exe.exists() and exe.stat().st_mtime > src.stat().st_mtime:
        return exe
    cmd = ["g++", "-O2", "-std=c++17", f"-I{ROOT}", str(src), "-lyaml-cpp", "-o", str(exe)]
    subprocess.run(cmd, check=True)
    return exe


def cantera_ref(mech: pathlib.Path, work: pathlib.Path) -> dict:
    """Cantera 参照を一時 CSV に出し、着火遅れ (max dT/dt) と終端温度だけを返す。"""
    import numpy as np, json
    csv = work / (mech.stem + "_ref.csv")
    subprocess.run([str(VENV_PY), str(ROOT / "tools/chem_reference_cantera.py"), str(mech), str(csv),
                    T0, P_ATM, T_END], check=True, capture_output=True)
    d = np.loadtxt(csv, delimiter=",", skiprows=1)
    t, T = d[:, 0], d[:, 1]
    i = int(np.argmax(np.diff(T) / np.maximum(np.diff(t), 1e-30)))
    return {"tau_ign_s": float(t[i + 1]), "T_end": float(T[-1])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--update", action="store_true", help="Cantera 参照を再生成する")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    GOLD.mkdir(exist_ok=True)
    import json
    gold_path = GOLD / "chem0d.json"
    gold = json.loads(gold_path.read_text()) if gold_path.exists() else {}
    work = pathlib.Path(os.environ.get("TMPDIR", "/tmp")) / "forge_chem0d"
    work.mkdir(parents=True, exist_ok=True)
    exe = build(work)
    ok = True
    for name, (mech, rtol) in MECHS.items():
        if a.update or name not in gold:
            if not VENV_PY.exists():
                print(f"[{name}] golden missing and no Cantera venv at {VENV_PY}"); ok = False; continue
            gold[name] = cantera_ref(mech, work)
            gold_path.write_text(json.dumps(gold, indent=1, sort_keys=True) + "\n")
            print(f"[{name}] golden regenerated: tau_ign {gold[name]['tau_ign_s']*1e6:.2f} us")
        r = subprocess.run([str(exe), str(mech), str(DB), "", T0, P_ATM, T_END],
                           env=ENV, capture_output=True, text=True)
        out = r.stdout + r.stderr
        jac = re.search(r"\[jacobian\].*\((PASS|FAIL)\)", out)
        tau = re.search(r"\[reactor abs-datum\].*T_end=([0-9.]+)K\s+tau_ign=([0-9.]+) us", out)
        mass = re.search(r"\[mass\]\s+Σω_s = ([0-9.e+-]+) kg/m3/s \(\|ω\|max ([0-9.e+-]+)\)", out)
        if r.returncode != 0 or not (jac and tau):
            print(f"[{name}] FAIL: test_chemistry did not produce results (rc={r.returncode})\n{out[-800:]}"); ok = False; continue
        tau_f = float(tau.group(2)) * 1e-6; tend_f = float(tau.group(1))
        tau_c = gold[name]["tau_ign_s"]; tend_c = gold[name]["T_end"]
        rel = abs(tau_f - tau_c) / tau_c
        mass_ok = (float(mass.group(1)) / float(mass.group(2)) < 1e-10) if mass else False
        verdict = "PASS" if (jac.group(1) == "PASS" and rel <= rtol and mass_ok and abs(tend_f - tend_c) < 10.0) else "FAIL"
        print(f"[{name}] {verdict}: tau_ign forge {tau_f*1e6:.2f} us vs Cantera {tau_c*1e6:.2f} us (rel {rel*100:.2f} %, tol {rtol*100:.0f} %), "
              f"T_end {tend_f:.0f} vs {tend_c:.0f} K, jacobian {jac.group(1)}, mass-conservation {'ok' if mass_ok else 'NG'}")
        ok &= (verdict == "PASS")
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
