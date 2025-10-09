#!/bin/bash

# Monitor memory usage of ZPGenHolo process
# Run this in one terminal while running test-large-zp.sh in another

echo "Waiting for ZPGenHolo to start..."
while ! pgrep -q ZPGenHolo; do
    sleep 0.5
done

sleep 0.5  # Give it a moment to fully start

# Get all PIDs, pick the one with highest memory (that's the actual binary)
PID=$(ps -eo pid,comm,rss | grep ZPGenHolo | grep -v grep | sort -k3 -rn | head -1 | awk '{print $1}')
if [ -z "$PID" ]; then
    echo "Error: Could not find ZPGenHolo process"
    exit 1
fi

echo "Found ZPGenHolo with PID: $PID"
echo "Monitoring memory usage (press Ctrl+C to stop)..."
echo ""
echo "Time(s)  | RSS(MB)  | VSZ(MB)"
echo "---------|----------|----------"

START_TIME=$(date +%s)
while ps -p $PID > /dev/null 2>&1; do
    CURRENT_TIME=$(date +%s)
    ELAPSED=$((CURRENT_TIME - START_TIME))

    # Get memory stats - ps returns in KB on macOS
    RSS_KB=$(ps -p $PID -o rss= 2>/dev/null | tr -d ' ')
    VSZ_KB=$(ps -p $PID -o vsz= 2>/dev/null | tr -d ' ')

    if [ ! -z "$RSS_KB" ]; then
        RSS_MB=$((RSS_KB / 1024))
        VSZ_MB=$((VSZ_KB / 1024))
        printf "%8d | %8d | %8d\n" $ELAPSED $RSS_MB $VSZ_MB
    fi

    sleep 1
done

echo ""
echo "Process finished"
