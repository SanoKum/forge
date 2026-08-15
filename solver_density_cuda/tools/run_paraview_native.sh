#!/usr/bin/env bash
# WSL/Linux native の ParaView を、forge の Python プラグインが読める状態で起動する。
#
# Ubuntu 24.04 の paraview 5.11.2 は Python 3.12 上で動くが、ParaView 同梱の
# paraview/detail/pythonalgorithm.py が Python 3.11 で削除された inspect.getargspec を
# import するため、Python プラグイン (tools/paraview/forge_filters.py) の読み込みが
#   WARN| Failed to load Python plugin: Failed to import `paraview.detail.pythonalgorithm`.
# で必ず失敗する。本スクリプトは tools/paraview/pvshim を PYTHONPATH に足して
# 起動時に getargspec を補うので、Tools > Manage Plugins からの読み込みも
# auto_load も通常どおり動く。
#
# 使い方:
#   ./tools/run_paraview_native.sh
#   ./tools/run_paraview_native.sh ../case/08.bump/run_0006_slau_loM_imp_cfl5/res_8000.xmf
#
# Docker 経由 (tools/run_paraview_gui.sh) の ParaView は Python 3.10 なので shim 不要。

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
shim_dir="$script_dir/paraview/pvshim"

if ! command -v paraview >/dev/null 2>&1; then
  echo "[run_paraview_native] ERROR: paraview が見つかりません (Docker 版は tools/run_paraview_gui.sh)" >&2
  exit 1
fi

export PYTHONPATH="$shim_dir${PYTHONPATH:+:$PYTHONPATH}"

args=()
for arg in "$@"; do
  # res_*.h5 を渡されたら、対応する xmf があればそちらを開く
  if [[ "$arg" == *.h5 && -f "${arg%.h5}.xmf" ]]; then
    args+=("${arg%.h5}.xmf")
  else
    args+=("$arg")
  fi
done

echo "[run_paraview_native] PYTHONPATH shim: $shim_dir" >&2
exec paraview "${args[@]}"
