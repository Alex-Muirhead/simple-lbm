echo "\n === CPU Performance ===\n"

echo "Runtime:"
for i in ./examples/taylor_green/*/cpu_performance.txt; do
    gawk -f ./tools/scrape_runtime.awk $i
done

echo "Speed:"
for i in ./examples/taylor_green/*/cpu_performance.txt; do
    gawk -f ./tools/scrape_speed.awk $i
done

echo "Bandwidth:"
for i in ./examples/taylor_green/*/cpu_performance.txt; do
    gawk -f ./tools/scrape_bandwidth.awk $i
done

echo "\n === GPU Performance ===\n"

echo "Runtime:"
for i in ./examples/taylor_green/*/gpu_performance.txt; do
    gawk -f ./tools/scrape_runtime.awk $i
done

echo "Speed:"
for i in ./examples/taylor_green/*/gpu_performance.txt; do
    gawk -f ./tools/scrape_speed.awk $i
done

echo "Bandwidth:"
for i in ./examples/taylor_green/*/gpu_performance.txt; do
    gawk -f ./tools/scrape_bandwidth.awk $i
done