# This isn't the Makefile for building the project. That's in Builds/LinuxMakefile. This is just utils.

build-decomp:
	cd build && cmake --build . --parallel

try: clean-logs
	cd build && ./decomp ../example-media/donut.wav ../example-media/talking.wav ../example-media/result.wav

debug: clean-logs
	cd build && gdb --args decomp ../example-media/donut.wav ../example-media/talking.wav ../example-media/result.wav

steps:
	./build/decomp --step HannWindow logs/020-carrier-highpass-b.txt step-logs/030-carrier-hann-b.txt

cross-steps:
	./build/decomp --step HannWindow logs/020-carrier-highpass-a.txt step-logs/030-carrier-hann-b.txt

clean-logs:
	rm logs/*-b.txt || echo "Logs clean already."

logs-to-sheet:
	node tools/logs-to-csv.js logs > logs/logs.csv

show-pd-snaps:
	node tools/wrap-images-in-html.js debug-snapshots/pd tools/page-template.html tools/entry-template.html > debug-snapshots/pd/index.html
	firefox debug-snapshots/pd/index.html

rsqrt: tools/try-rsqrt
	./tools/try-rsqrt logs/050-carrier-rfft-added-a.txt > logs/060-carrier-rsqrt-experimental.txt

tools/try-rsqrt:
	g++ tools/try-rsqrt.c -o tools/try-rsqrt -std=c++11

build-try-rsqrt:
	g++ tools/try-rsqrt.c -o tools/try-rsqrt -std=c++11
