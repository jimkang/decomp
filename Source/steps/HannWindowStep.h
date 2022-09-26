//
// Created by jimkang on 3/16/22.
//

#ifndef DECOMP_HANNWINDOWSTEP_H
#define DECOMP_HANNWINDOWSTEP_H

#include "JuceHeader.h"
#include <math.h>
#include "../Consts.h"

static void applyHannWindow(float *signalBlock, int size) {
    juce::dsp::WindowingFunction<float> window(size, juce::dsp::WindowingFunction<float>::hann, false);
    window.multiplyWithWindowingTable(signalBlock, size);

    //float maxOrig = 0.0;
    //float maxNew = 0.0;
//
    //// TODO: Hann curve should be memoized.
    //for (int i = 0; i < size; ++i) {
    //const float orig = signalBlock[i];
    //const float hann = (cos( (i*1.0/size * 2.0 + 1.0)*M_PI ) + 1)/2.0;
    //signalBlock[i] = orig * hann;
    //if (orig > maxOrig) {
    //maxOrig = orig;
    //}
    //if (signalBlock[i] > maxNew) {
    //maxNew = signalBlock[i];
    //}
    //}
//
    //cout << "maxOrig: " << maxOrig << ", maxNew: " << maxNew << endl;
}


#endif //DECOMP_HANNWINDOWSTEP_H
