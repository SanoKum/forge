#!/usr/bin/env python3
"""PreToolUse(Bash) フック: forge の直接実行を弾き run_case.sh 経由を強制する。

stdin に Claude Code のフック入力 JSON を受け、`build/forge` を実行しようとする
Bash コマンド (run_case.sh 経由でない、かつ freshness 用 find -newer でない) を deny する。
deny は permissionDecision で返す。該当しなければ {} (許可)。

注意: シェル変数展開 ($FORGE 等) で隠れた呼び出しは検知できない。最終的な担保は
Stop フック (hook_convergence_gate.py, 実行方法に依らずファイル状態で検査) が行う。
"""
import sys, json, re

try:
    data = json.load(sys.stdin)
except Exception:
    print("{}"); sys.exit(0)

cmd = ((data.get("tool_input") or {}).get("command") or "")

exec_forge = re.search(r'(^|[;&|]|\bnohup\s+|bash\s+-c\s+["\']?)\s*\S*build/forge\b', cmd) is not None
benign = ("run_case.sh" in cmd) or ("-newer" in cmd) or (re.search(r'\bfind\b', cmd) is not None)

if exec_forge and not benign:
    out = {"hookSpecificOutput": {
        "hookEventName": "PreToolUse",
        "permissionDecision": "deny",
        "permissionDecisionReason": (
            "forge は solver_density_cuda/tools/run_case.sh 経由で実行すること。"
            "ラッパーが実行後に check_convergence.py を走らせ CONVERGENCE_VERDICT.txt を残すため、"
            "収束確認が自動で担保される。直接 build/forge を呼ばない "
            "(例: `nohup solver_density_cuda/tools/run_case.sh <run_dir> &`)。")}}
    print(json.dumps(out)); sys.exit(0)

print("{}")
