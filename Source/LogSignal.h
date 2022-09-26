#include <iostream>
#include <fstream>
#include "Consts.h"

void logSignalToPath(int signalLength, const float *signal, const string &filePath);

#pragma once

using namespace std;

// TODO: Clear logs at beginning of program.
void logSignal(const char *fileName, int signalLength, const float *signal) {
    string filePath = baseLogPath;
    filePath += fileName;
    logSignalToPath(signalLength, signal, filePath);
}

void logSignalToPath(int signalLength, const float *signal, const string &filePath) {
    ofstream out(filePath.c_str(), ios_base::app);
    out.precision(14);
    for (int i = 0; i < signalLength; ++i) {
        out << signal[i] << endl;
    }
    out.close();
}
