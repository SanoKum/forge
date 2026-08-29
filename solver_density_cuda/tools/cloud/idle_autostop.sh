#!/usr/bin/env bash
# GPU アイドルが続いたらインスタンスを自動 stop する消し忘れ保険。
# root の cron から 5 分おきに呼ばれる想定 (setup_instance.sh が /etc/cron.d に登録)。
# 「GPU 使用率 0 かつ forge プロセス無し かつ ログインセッション無し」が
# IDLE_LIMIT 回 (既定 6 回 = 30 分) 連続したら shutdown する。
set -u

IDLE_LIMIT="${FORGE_IDLE_LIMIT:-6}"
state_file="/var/run/forge-idle-count"

gpu_util="$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null | sort -rn | head -1)"
forge_running="$(pgrep -c -x forge || true)"
sessions="$(who | wc -l)"

if [[ -z "$gpu_util" ]]; then
  # nvidia-smi 失敗時は判定不能として何もしない
  exit 0
fi

if [[ "$gpu_util" -eq 0 && "$forge_running" -eq 0 && "$sessions" -eq 0 ]]; then
  count=$(( $(cat "$state_file" 2>/dev/null || echo 0) + 1 ))
else
  count=0
fi
echo "$count" > "$state_file"

if [[ "$count" -ge "$IDLE_LIMIT" ]]; then
  echo "$(date -Is) idle x${count} -> shutdown"
  shutdown -h now
fi
