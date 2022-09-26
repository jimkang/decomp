#include <JuceHeader.h>
#include <math.h>
#include "Utils.h"
#include "Consts.h"
#include "LogSignal.h"
#include "Sqrt.h"
#include "steps/HannWindowStep.h"

using namespace juce;
using namespace std;

static void
decompChannel(vector<float> &carrierSamples, vector<float> &infoSamples, double sampleRate, vector<float> &outSamples);

static void
decompBlock(vector<float> &carrierBlockSamples, vector<float> &infoBlockSamples, IIRFilter &carrierHighPassFilter,
            IIRFilter &infoHighPassFilter, vector<float> &outBlockSamples, int blockIndex);

static void getReducedCombinedAmpFactors(
        ComplexFFTArray &carrierFFTData, ComplexFFTArray &infoFFTData, FFTArray &reducedAmpFactors);

static void decomp(AudioBuffer<float> &carrierBuffer, AudioBuffer<float> &infoBuffer,
                   double sampleRate,
                   AudioBuffer<float> &outBuffer) {

    int channelCount = carrierBuffer.getNumChannels();
    int infoBufferChannelCount = infoBuffer.getNumChannels();
    if (infoBufferChannelCount < channelCount) {
        channelCount = infoBufferChannelCount;
    }
    for (int ch = 0; ch < channelCount; ++ch) {
        const int maxSamples = outBuffer.getNumSamples();
        // Copy the arrays into vectors.
        vector<float> carrierSamples(carrierBuffer.getReadPointer(ch), carrierBuffer.getReadPointer(ch) + maxSamples);
        vector<float> infoSamples(infoBuffer.getReadPointer(ch), infoBuffer.getReadPointer(ch) + maxSamples);
        vector<float> outSamples(maxSamples);

        //for (int i = 0; i < outSamples.size(); ++i) {
        //outSamples[i] = carrierSamples[i];
        //}

        decompChannel(carrierSamples, infoSamples, sampleRate, outSamples);

        //for (auto it = outSamples.begin(); it != outSamples.end(); ++it) {
        //if (*it > 0.01) {
        //cout << "hey";
        //}
        //}
        float *outWritePtr = outBuffer.getWritePointer(ch);
        //cout << "outSamples size: " << outSamples.size() << endl;
        for (int i = 0; i < outSamples.size(); ++i) {
            outWritePtr[i] = outSamples[i];
            if (outSamples[i] > 0) {
                int x = 2;
            }
        }
    }
}

static void
decompChannel(vector<float> &carrierSamples, vector<float> &infoSamples, double sampleRate, vector<float> &outSamples) {
    const int maxBlocks = outSamples.size() / blockSize;
    // Leave out the last partial block for now.

    auto carrierStart = carrierSamples.begin();
    auto infoStart = infoSamples.begin();
    auto outStart = outSamples.begin();

    IIRFilter carrierHighPassFilter;
    IIRFilter infoHighPassFilter;
    double freq = 5;
    carrierHighPassFilter.setCoefficients(IIRCoefficients::makeHighPass(sampleRate, freq));
    infoHighPassFilter.setCoefficients(IIRCoefficients::makeHighPass(sampleRate, freq));

    //for (int i = 0; i < outSamples.size(); ++i) {
    //outSamples[i] = carrierSamples[i];
    //}
    //return;

    for (int blockIndex = 0; blockIndex < maxBlocks; ++blockIndex) {
        cout << "vocoding block  " << blockIndex << endl;
        const int sampleIndex = blockIndex * blockSize;

        vector<float> carrierBlockSamples(blockSize);
        vector<float> infoBlockSamples(blockSize);
        vector<float> outBlockSamples(blockSize);

        buildBlockWithOverlap(carrierSamples.data(), sampleIndex, blockSize, overlapFactor, carrierBlockSamples);
        buildBlockWithOverlap(infoSamples.data(), sampleIndex, blockSize, overlapFactor, infoBlockSamples);

        decompBlock(
                carrierBlockSamples,
                infoBlockSamples,
                carrierHighPassFilter,
                infoHighPassFilter,
                outBlockSamples,
                blockIndex
        );

        // TODO: Cut down on the copying.
        for (int outBlockSampleIndex = 0; outBlockSampleIndex < blockSize; ++outBlockSampleIndex) {
            *outStart = outBlockSamples[outBlockSampleIndex];
            ++outStart;
        }
    }
}

static void
decompBlock(vector<float> &carrierBlockSamples, vector<float> &infoBlockSamples, IIRFilter &carrierHighPassFilter,
            IIRFilter &infoHighPassFilter, vector<float> &outBlockSamples, int blockIndex) {

    auto carrierSample = carrierBlockSamples.begin();
    //auto outSample = outBlockSamples.begin();
    //for (int i = 0; i < outBlockSamples.size(); ++i) {
    //outBlockSamples[i] = carrierSample[i];
    //}
    //return;

    if (blockIndex == blockIndexToLog) {
        logSignal("005-info-raw-b.txt", infoBlockSamples.size(), infoBlockSamples.data());
        logSignal("010-carrier-raw-b.txt", carrierBlockSamples.size(), carrierBlockSamples.data());
    }

    // Drop high-pass for now.
    //carrierHighPassFilter.processSamples(carrierBlockSamples.data(), carrierBlockSamples.size());
    //infoHighPassFilter.processSamples(infoBlockSamples.data(), infoBlockSamples.size());
    if (blockIndex == blockIndexToLog) {
        logSignal("006-info-highpass-b.txt", infoBlockSamples.size(), infoBlockSamples.data());
        logSignal("020-carrier-highpass-b.txt", carrierBlockSamples.size(), carrierBlockSamples.data());
    }

    applyHannWindow(infoBlockSamples.data(), infoBlockSamples.size());
    applyHannWindow(carrierBlockSamples.data(), carrierBlockSamples.size());

    // TODO: Include channel in filename.
    if (blockIndex == blockIndexToLog) {
        logSignal("007-info-hann-b.txt", infoBlockSamples.size(), infoBlockSamples.data());
        logSignal("030-carrier-hann-b.txt", carrierBlockSamples.size(), carrierBlockSamples.data());
    }

    // Run a real-only FFT on both signals.
    ComplexFFTArray carrierFFTData;
    ComplexFFTArray infoFFTData;

    for (int i = 0; i < carrierBlockSamples.size(); ++i) {
        carrierFFTData[i] = carrierBlockSamples[i];
        infoFFTData[i] = infoBlockSamples[i];
    }

    getFFT(carrierFFTData);
    getFFT(infoFFTData);

    //if (blockIndex == blockIndexToLog) {
    //logSignal("040-carrierFFT-b.txt", fftSize, carrierFFTData.data());
    //}

    FFTArray infoRealBins;
    getReal(infoFFTData, infoRealBins);
    if (blockIndex == blockIndexToLog) {
        logSignal("042-info-fft-real-b.txt", fftSize, infoRealBins.data());
    }


    FFTArray infoImagBins;
    getImaginary(infoFFTData, infoImagBins);
    if (blockIndex == blockIndexToLog) {
        logSignal("043-info-fft-imag-b.txt", fftSize, infoImagBins.data());
    }

    // Multiply the reduced real components of the carrier fft by the reduced
    // combined amps.
    FFTArray carrierRealBins;
    getReal(carrierFFTData, carrierRealBins);
    if (blockIndex == blockIndexToLog) {
        logSignal("045-carrier-fft-real-b.txt", fftSize, carrierRealBins.data());
    }

    FFTArray carrierImagBins;
    getImaginary(carrierFFTData, carrierImagBins);
    if (blockIndex == blockIndexToLog) {
        logSignal("047-carrier-fft-imag-b.txt", fftSize, carrierImagBins.data());
    }

    // NOTE: we are altering the things we square.
    squareSignal(carrierRealBins.data(), fftSize);
    squareSignal(carrierImagBins.data(), fftSize);
    squareSignal(infoRealBins.data(), fftSize);
    squareSignal(infoImagBins.data(), fftSize);
    if (blockIndex == blockIndexToLog) {
        logSignal("048-info-fft-real-sq-b.txt", fftSize, infoRealBins.data());
        logSignal("048.5-info-fft-imag-sq-b.txt", fftSize, infoImagBins.data());
        logSignal("049-carrier-fft-real-sq-b.txt", fftSize, carrierRealBins.data());
        logSignal("049.5-carrier-fft-imag-sq-b.txt", fftSize, carrierImagBins.data());
    }

    FFTArray carrierFFTSqAdded;
    FFTArray infoFFTSqAdded;
    FloatVectorOperations::add(infoFFTSqAdded.data(), infoRealBins.data(), infoImagBins.data(), fftSize);
    FloatVectorOperations::add(carrierFFTSqAdded.data(), carrierRealBins.data(), carrierImagBins.data(), fftSize);
    if (blockIndex == blockIndexToLog) {
        logSignal("050-carrier-rfft-added-b.txt", fftSize, carrierFFTSqAdded.data());
        logSignal("055-info-rfft-added-b.txt", fftSize, infoFFTSqAdded.data());
    }

    FFTArray carrierFFTSqAddedRSqrt;
    FFTArray infoFFTSqAddedSqrt;
    // Why does this get so big?
    rSqrtSignal(carrierFFTSqAdded.data(), fftSize, carrierFFTSqAddedRSqrt.data());
    sqrtSignal(infoFFTSqAdded.data(), fftSize, infoFFTSqAddedSqrt.data());
    if (blockIndex == blockIndexToLog) {
        logSignal("060-carrier-rsqrt-b.txt", fftSize, carrierFFTSqAddedRSqrt.data());
        logSignal("070-info-sqrt-b.txt", fftSize, infoFFTSqAddedSqrt.data());
    }

    FFTArray combinedAmpFactors;
    FloatVectorOperations::multiply(
            combinedAmpFactors.data(), // dest
            carrierFFTSqAddedRSqrt.data(),
            infoFFTSqAddedSqrt.data(),
            fftSize);
    //for (int i = 0; i < combinedAmpFactors.size(); ++i) {
    //combinedAmpFactors[i] = carrierFFTSqAddedRSqrt[i] * infoFFTSqAddedSqrt[i];
    //}
    if (blockIndex == blockIndexToLog) {
        logSignal("080-amp-factor-roots-multiplied-b.txt", fftSize, combinedAmpFactors.data());
    }

    // Turn down the combined amps.
    FFTArray reducedAmpFactors;
    FloatVectorOperations::multiply(
            reducedAmpFactors.data(),
            combinedAmpFactors.data(),
            1.0 / hannOverlapGain,// * smallifyFactor,
            fftSize);

    FFTArray carrierRealWithReducedAmpFactors;
    FloatVectorOperations::multiply(
            carrierRealWithReducedAmpFactors.data(),
            carrierRealBins.data(),
            reducedAmpFactors.data(),
            fftSize);
    if (blockIndex == blockIndexToLog) {
        logSignal("100-carrier-fft-real-x-reduced-amp-factors-b.txt", fftSize, carrierRealWithReducedAmpFactors.data());
    }

    // Multiply the imaginary components of the carrier fft by the reduced
    // combined amps.
    FFTArray carrierImagWithReducedAmpFactors;
    FloatVectorOperations::multiply(
            carrierImagWithReducedAmpFactors.data(),
            carrierImagBins.data(),
            reducedAmpFactors.data(),
            fftSize);
    if (blockIndex == blockIndexToLog) {
        logSignal("150-carrier-fft-imag-x-reduced-amp-factors-b.txt", fftSize, carrierImagWithReducedAmpFactors.data());
    }

    ComplexFFTArray ifftData;
    getIFFT(carrierRealWithReducedAmpFactors, carrierImagWithReducedAmpFactors, ifftData);
    //getIFFT(carrierRealBins, carrierImagBins, ifftData);

    applyHannWindow(ifftData.data(), ifftData.size());

    // Copy the results to the channel.
    for (int i = 0; i < fftSize; ++i) {
        outBlockSamples[i] = ifftData[i];
    }
}
