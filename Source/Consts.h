#pragma once

#include <math.h>

const int fftPowerOf2 = 10;
const int fftSize = 1 << fftPowerOf2;
const int blockSize = fftSize;
const int blockSizeInBytes = blockSize * sizeof(float);
// These come from Miller Puckette's timbre-stamp Pure Data patch.
// I imagine they were determined via experimentation.
const float maxCarrierMag = 9;
const float smallifyFactor = 0.00065;
const float tinyNumber = pow(10, -20);

const int overlapFactor = 4;
const int overlapOffset = blockSize / overlapFactor;
// Every time the Hann window is applied, the amplitude is reduced to
// half. (We are applying it twice.)
// However, we are also adding the signal to itself (albeit with
// offsets), overlapFactor number of times.
const float hannOverlapGain = pow(0.5, 2) * overlapFactor;

// This is to help approximate what Pure Data's sqrt does.
const float closeEnoughToZero = pow(10, -8.5);

const float minimumNonZeroFFTResult = 3E-05;

// TODO: Make this an arg.
const char *baseLogPath = "../logs/";
const int blockIndexToLog = 68;
