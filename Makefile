# This isn't the Makefile for building the project. That's in Builds/LinuxMakefile. This is just utils.

build-decomp:
	cd build && cmake --build . --parallel

try:
	cd build && ./decomp ../example-media/donut.wav 8 ../example-media

debug:
	cd build && gdb --args decomp ../example-media/donut.wav 8 ../example-media/

steps:
	./build/decomp --step HannWindow logs/020-carrier-highpass-b.txt step-logs/030-carrier-hann-b.txt

cross-steps:
	./build/decomp --step HannWindow logs/020-carrier-highpass-a.txt step-logs/030-carrier-hann-b.txt
