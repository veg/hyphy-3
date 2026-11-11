#!/usr/bin/env bash
# Manage the zero-account WALRUS Cluster Gateway + Cloudflared HTTPS Tunnel on silverback.
# Usage: ./manage_service.sh {start|stop|restart|status|url} [named_tunnel_token_optional]
set -euo pipefail

SERVICE_ROOT="${WALRUS_SERVICE_ROOT:-$HOME/walrus-service}"
GATEWAY_DIR="$SERVICE_ROOT/gateway"
LOG_DIR="$SERVICE_ROOT/logs"
CLOUDFLARED_BIN="$SERVICE_ROOT/bin/cloudflared"
BIND_PORT="${WALRUS_BIND_PORT:-18791}"
SESSION="walrus_service"

mkdir -p "$GATEWAY_DIR" "$LOG_DIR" "$SERVICE_ROOT/jobs"

cmd="${1:-status}"
named_token="${2:-${CLOUDFLARED_TUNNEL_TOKEN:-}}"

extract_url() {
  if [ -f "$LOG_DIR/cloudflared.log" ]; then
    grep -Eo 'https://[a-z0-9-]+\.trycloudflare\.com' "$LOG_DIR/cloudflared.log" | tail -n 1 || true
  fi
}

case "$cmd" in
  start)
    if tmux has-session -t "$SESSION" 2>/dev/null; then
      echo "Session '$SESSION' is already running."
    else
      rm -f "$LOG_DIR/cloudflared.log"
      tmux new-session -d -s "$SESSION" -n gateway \
        "bash -lc 'while true; do /usr/bin/python3 $GATEWAY_DIR/walrus_gateway.py >> $LOG_DIR/gateway.log 2>&1; sleep 2; done'"

      if [ -n "$named_token" ]; then
        tmux new-window -t "$SESSION" -n tunnel \
          "bash -lc 'while true; do $CLOUDFLARED_BIN tunnel --no-autoupdate run --token \"$named_token\" >> $LOG_DIR/cloudflared.log 2>&1; sleep 3; done'"
      else
        tmux new-window -t "$SESSION" -n tunnel \
          "bash -lc 'while true; do $CLOUDFLARED_BIN tunnel --no-autoupdate --url http://127.0.0.1:$BIND_PORT >> $LOG_DIR/cloudflared.log 2>&1; sleep 3; done'"
      fi
      echo "Started tmux session '$SESSION' (gateway + cloudflared)."
    fi

    # Wait up to 12s for quick-tunnel URL if applicable
    if [ -z "$named_token" ]; then
      for _ in $(seq 1 24); do
        url="$(extract_url)"
        if [ -n "$url" ]; then
          echo "$url" > "$LOG_DIR/tunnel_url.txt"
          echo "Live HTTPS Endpoint: $url"
          exit 0
        fi
        sleep 0.5
      done
    fi
    ;;

  stop)
    if tmux has-session -t "$SESSION" 2>/dev/null; then
      tmux kill-session -t "$SESSION"
      echo "Stopped '$SESSION'."
    else
      echo "Session '$SESSION' is not running."
    fi
    ;;

  restart)
    "$0" stop
    sleep 1
    "$0" start "$named_token"
    ;;

  url)
    url="$(extract_url)"
    if [ -n "$url" ]; then
      echo "$url"
    elif [ -f "$LOG_DIR/tunnel_url.txt" ]; then
      cat "$LOG_DIR/tunnel_url.txt"
    else
      echo "No Quick Tunnel URL detected (check $LOG_DIR/cloudflared.log or named tunnel DNS)." >&2
      exit 1
    fi
    ;;

  status)
    if tmux has-session -t "$SESSION" 2>/dev/null; then
      echo "Service tmux session '$SESSION': RUNNING"
      tmux list-windows -t "$SESSION"
      curl -fsS "http://127.0.0.1:$BIND_PORT/v1/status" && echo "" || echo "Local HTTP gateway not responding yet."
      url="$(extract_url)"
      [ -n "$url" ] && echo "Public Quick Tunnel URL: $url"
    else
      echo "Service tmux session '$SESSION': STOPPED"
      exit 1
    fi
    ;;

  *)
    echo "Usage: $0 {start|stop|restart|status|url} [named_tunnel_token]" >&2
    exit 2
    ;;
esac
