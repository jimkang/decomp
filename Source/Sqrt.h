#include <math.h>
#include "Consts.h"

static float reciprocalSqRt(float bin) {
    return 1.0 / sqrt(bin);
}

static void rSqrtSignal(const float *array, int size, float *outArray) {
    for (int i = 0; i < size; ++i) {
        float val = array[i];
        val += tinyNumber;
        const float result = reciprocalSqRt(val);
        outArray[i] = result;
    }
}

static void sqrtSignal(const float *array, int size, float *outArray) {
    for (int i = 0; i < size; ++i) {
        // TODO: Move this out of here and just add
        // it to the incoming signal ahead of time.
        float val = array[i];
        if (val == 0) {
            val = tinyNumber;
        }
        outArray[i] = sqrt(val);
    }
}
