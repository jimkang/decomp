#ifndef DECOMP_RUNSTEP_H
#define DECOMP_RUNSTEP_H

#endif //DECOMP_RUNSTEP_H

#include <ios>
#include <iostream>
#include <fstream>
#include <vector>
#include "steps/HannWindowStep.h"

using namespace std;

int runStep(const char *stepName, const char *inTextPath, const char *outTextPath) {
    ifstream in(inTextPath, ios_base::app);
    if (!in.is_open()) {
        std::cout << "failed to open " << inTextPath << endl;
        return 1;
    }

    in.precision(14);
    vector<float> values;
    string line;
    while (getline(in, line)) {
        float value = strtof(line.c_str(), nullptr);
        values.push_back(value);
    }
    in.close();

//    for (int i = 0; i < values.size(); ++i) {
//        cout << "value:" << values[i] << endl;
//    }

    if (strcmp(stepName, "HannWindow") == 0) {
        applyHannWindow(values.data(), (int)values.size());
    }

    logSignalToPath((int)values.size(), values.data(), outTextPath);

    return 0;
}