#!/usr/bin/env python3
"""Stop フック: forge を回したのに収束チェックしていない run があればターン終了を block する。

case/*/run_*/residual_history.csv のうち「最近 (既定 180 分) 更新された」ものを対象に、
同じ run の CONVERGENCE_VERDICT.txt が無い / residual より古い ものを検出したら block する。
これは forge の実行方法 (直接/ラッパー/バックグラウンド) に依らずファイル状態だけで判定するので、
「位置安定だけ見て収束と判断する」近道を構造的に防ぐ。

block 時は decision=block + reason を返し、解消法 (run_case.sh 経由 or check_convergence.py で
VERDICT を残す) と「収束/一致を主張するなら VERDICT 行を引用」を伝える。
"""
import sys, os, glob, json, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
RECENT_SEC = 180 * 60

bad = []
now = time.time()
for rh in glob.glob(os.path.join(ROOT, "case", "*", "run_*", "residual_history.csv")):
    try:
        m = os.path.getmtime(rh)
    except OSError:
        continue
    if now - m > RECENT_SEC:
        continue  # 旧 run は対象外 (このセッションの活動のみ強制)
    v = os.path.join(os.path.dirname(rh), "CONVERGENCE_VERDICT.txt")
    if (not os.path.exists(v)) or (os.path.getmtime(v) < m):
        bad.append(os.path.relpath(os.path.dirname(rh), ROOT))

if bad:
    reason = (
        "forge を実行したのに収束チェック未実施/古い run があります:\n  - "
        + "\n  - ".join(sorted(bad))
        + "\n\n各 run を `solver_density_cuda/tools/run_case.sh <run_dir>` 経由で実行する"
        "(自動で VERDICT を残す)か、`python3 solver_density_cuda/tools/check_convergence.py "
        "<run_dir> > <run_dir>/CONVERGENCE_VERDICT.txt` を実行して VERDICT を残すこと。\n"
        "**収束/一致/ショックレス等を主張する応答では check_convergence.py の VERDICT 行を必ず引用する。"
        "衝撃波位置の安定や rms_ro 単独で『収束』と判断しない。**")
    print(json.dumps({"decision": "block", "reason": reason}))
else:
    print("{}")
