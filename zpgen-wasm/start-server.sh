#!/bin/bash
# Start web server for ZPGen WASM browser interface

PORT=${1:-8888}
DIR="examples/browser"

echo "Starting ZPGen WASM web server..."
echo "Directory: $DIR"
echo "Port: $PORT"
echo ""
echo "Open in browser: http://localhost:$PORT"
echo ""
echo "Press Ctrl+C to stop"
echo ""

cd "$DIR" && python3 -m http.server "$PORT"
