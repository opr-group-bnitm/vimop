#!/usr/bin/env bash
set -u
pid="${1:?pid required}"
log="${2:?log required}"
while kill -0 "$pid" 2>/dev/null; do
  {
    printf '\n[%s]\n' "$(date -Is)"
    ps -p "$pid" -o pid=,etime=,pcpu=,pmem=,stat=,cmd=
    df -h /data /gemini
  } >> "$log"
  sleep 60
done
printf '\n[%s] process exited\n' "$(date -Is)" >> "$log"
