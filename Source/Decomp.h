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
decompChannel(vector<float> &inputSamples, vector<float> &infoSamples, double sampleRate, vector<float> &outSamples);

static void
decompBlock(vector<float> &inputBlockSamples, vector<float> &infoBlockSamples, IIRFilter &inputHighPassFilter,
            IIRFilter &infoHighPassFilter, vector<float> &outBlockSamples, int blockIndex);

static void getReducedCombinedAmpFactors(
        ComplexFFTArray &inputFFTData, ComplexFFTArray &infoFFTData, FFTArray &reducedAmpFactors);

static void decomp(AudioBuffer<float> &inputBuffer,
                   double sampleRate,
                   int bucketCount,
                   vector<juce::AudioBuffer<float>> &outBuffers) {

  int channelCount = inputBuffer.getNumChannels();
  int inputBufferChannelCount = inputBuffer.getNumChannels();
  if (inputBufferChannelCount < channelCount) {
      channelCount = inputBufferChannelCount;
  }
  const int maxSamples = inputBuffer.getNumSamples();
  for (int ch = 0; ch < channelCount; ++ch) {
    // Copy the arrays into vectors.
    vector<float> inputSamples(inputBuffer.getReadPointer(ch), inputBuffer.getReadPointer(ch) + maxSamples);
    vector<vector<float>> outSamplesCollection;

    for (auto it = outBuffers.begin(); it != outBuffers.end(); ++it) {
      vector<float> outputSamples(
        it->getReadPointer(ch), it->getReadPointer(ch) + maxSamples
      );

      // Test by copying.
      for (int i = 0; i < outputSamples.size(); ++i) {
        outputSamples[i] = inputSamples[i];
      }
      cout << "outputSamples value: " << outputSamples[1000] << endl;

      // Whoa, this is copying the vector. If you do this before setting the
      // values, the changes don't make it to outSamplesCollection!
      outSamplesCollection.push_back(outputSamples);
    }

    //decompChannel(inputSamples, infoSamples, sampleRate, outSamples);

    //for (auto it = outSamples.begin(); it != outSamples.end(); ++it) {
    //if (*it > 0.01) {
    //cout << "hey";
    //}
    //}
    int bufferNumber = 0;
    for (auto it = outBuffers.begin(); it != outBuffers.end(); ++it) {
      float *outWritePtr = it->getWritePointer(ch);
      vector<float> outputSamples = outSamplesCollection[bufferNumber];
      cout << "outputSamples value at end: " << outputSamples[1000] << endl;
      for (int i = 0; i < outputSamples.size(); ++i) {
        outWritePtr[i] = outputSamples[i];
      }
      bufferNumber += 1;
    }
  }
}

static void
decompChannel(vector<float> &inputSamples, vector<float> &infoSamples, double sampleRate, vector<float> &outSamples) {
    const int maxBlocks = outSamples.size() / blockSize;
    // Leave out the last partial block for now.

    auto inputStart = inputSamples.begin();
    auto infoStart = infoSamples.begin();
    auto outStart = outSamples.begin();

    IIRFilter inputHighPassFilter;
    IIRFilter infoHighPassFilter;
    double freq = 5;
    inputHighPassFilter.setCoefficients(IIRCoefficients::makeHighPass(sampleRate, freq));
    infoHighPassFilter.setCoefficients(IIRCoefficients::makeHighPass(sampleRate, freq));

    //for (int i = 0; i < outSamples.size(); ++i) {
    //outSamples[i] = inputSamples[i];
    //}
    //return;

    for (int blockIndex = 0; blockIndex < maxBlocks; ++blockIndex) {
        cout << "vocoding block  " << blockIndex << endl;
        const int sampleIndex = blockIndex * blockSize;

        vector<float> inputBlockSamples(blockSize);
        vector<float> infoBlockSamples(blockSize);
        vector<float> outBlockSamples(blockSize);

        buildBlockWithOverlap(inputSamples.data(), sampleIndex, blockSize, overlapFactor, inputBlockSamples);
        buildBlockWithOverlap(infoSamples.data(), sampleIndex, blockSize, overlapFactor, infoBlockSamples);

        decompBlock(
                inputBlockSamples,
                infoBlockSamples,
                inputHighPassFilter,
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
decompBlock(vector<float> &inputBlockSamples, vector<float> &infoBlockSamples, IIRFilter &inputHighPassFilter,
            IIRFilter &infoHighPassFilter, vector<float> &outBlockSamples, int blockIndex) {

    auto inputSample = inputBlockSamples.begin();
    //auto outSample = outBlockSamples.begin();
    //for (int i = 0; i < outBlockSamples.size(); ++i) {
    //outBlockSamples[i] = inputSample[i];
    //}
    //return;

    if (blockIndex == blockIndexToLog) {
        logSignal("005-info-raw-b.txt", infoBlockSamples.size(), infoBlockSamples.data());
        logSignal("010-input-raw-b.txt", inputBlockSamples.size(), inputBlockSamples.data());
    }

    // Drop high-pass for now.
    //inputHighPassFilter.processSamples(inputBlockSamples.data(), inputBlockSamples.size());
    //infoHighPassFilter.processSamples(infoBlockSamples.data(), infoBlockSamples.size());
    if (blockIndex == blockIndexToLog) {
        logSignal("006-info-highpass-b.txt", infoBlockSamples.size(), infoBlockSamples.data());
        logSignal("020-input-highpass-b.txt", inputBlockSamples.size(), inputBlockSamples.data());
    }

    applyHannWindow(infoBlockSamples.data(), infoBlockSamples.size());
    applyHannWindow(inputBlockSamples.data(), inputBlockSamples.size());

    // TODO: Include channel in filename.
    if (blockIndex == blockIndexToLog) {
        logSignal("007-info-hann-b.txt", infoBlockSamples.size(), infoBlockSamples.data());
        logSignal("030-input-hann-b.txt", inputBlockSamples.size(), inputBlockSamples.data());
    }

    // Run a real-only FFT on both signals.
    ComplexFFTArray inputFFTData;
    ComplexFFTArray infoFFTData;

    for (int i = 0; i < inputBlockSamples.size(); ++i) {
        inputFFTData[i] = inputBlockSamples[i];
        infoFFTData[i] = infoBlockSamples[i];
    }

    getFFT(inputFFTData);
    getFFT(infoFFTData);

    //if (blockIndex == blockIndexToLog) {
    //logSignal("040-inputFFT-b.txt", fftSize, inputFFTData.data());
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

    // Multiply the reduced real components of the input fft by the reduced
    // combined amps.
    FFTArray inputRealBins;
    getReal(inputFFTData, inputRealBins);
    if (blockIndex == blockIndexToLog) {
        logSignal("045-input-fft-real-b.txt", fftSize, inputRealBins.data());
    }

    FFTArray inputImagBins;
    getImaginary(inputFFTData, inputImagBins);
    if (blockIndex == blockIndexToLog) {
        logSignal("047-input-fft-imag-b.txt", fftSize, inputImagBins.data());
    }

    // NOTE: we are altering the things we square.
    squareSignal(inputRealBins.data(), fftSize);
    squareSignal(inputImagBins.data(), fftSize);
    squareSignal(infoRealBins.data(), fftSize);
    squareSignal(infoImagBins.data(), fftSize);
    if (blockIndex == blockIndexToLog) {
        logSignal("048-info-fft-real-sq-b.txt", fftSize, infoRealBins.data());
        logSignal("048.5-info-fft-imag-sq-b.txt", fftSize, infoImagBins.data());
        logSignal("049-input-fft-real-sq-b.txt", fftSize, inputRealBins.data());
        logSignal("049.5-input-fft-imag-sq-b.txt", fftSize, inputImagBins.data());
    }

    FFTArray inputFFTSqAdded;
    FFTArray infoFFTSqAdded;
    FloatVectorOperations::add(infoFFTSqAdded.data(), infoRealBins.data(), infoImagBins.data(), fftSize);
    FloatVectorOperations::add(inputFFTSqAdded.data(), inputRealBins.data(), inputImagBins.data(), fftSize);
    if (blockIndex == blockIndexToLog) {
        logSignal("050-input-rfft-added-b.txt", fftSize, inputFFTSqAdded.data());
        logSignal("055-info-rfft-added-b.txt", fftSize, infoFFTSqAdded.data());
    }

    FFTArray inputFFTSqAddedRSqrt;
    FFTArray infoFFTSqAddedSqrt;
    // Why does this get so big?
    rSqrtSignal(inputFFTSqAdded.data(), fftSize, inputFFTSqAddedRSqrt.data());
    sqrtSignal(infoFFTSqAdded.data(), fftSize, infoFFTSqAddedSqrt.data());
    if (blockIndex == blockIndexToLog) {
        logSignal("060-input-rsqrt-b.txt", fftSize, inputFFTSqAddedRSqrt.data());
        logSignal("070-info-sqrt-b.txt", fftSize, infoFFTSqAddedSqrt.data());
    }

    FFTArray combinedAmpFactors;
    FloatVectorOperations::multiply(
            combinedAmpFactors.data(), // dest
            inputFFTSqAddedRSqrt.data(),
            infoFFTSqAddedSqrt.data(),
            fftSize);
    //for (int i = 0; i < combinedAmpFactors.size(); ++i) {
    //combinedAmpFactors[i] = inputFFTSqAddedRSqrt[i] * infoFFTSqAddedSqrt[i];
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

    FFTArray inputRealWithReducedAmpFactors;
    FloatVectorOperations::multiply(
            inputRealWithReducedAmpFactors.data(),
            inputRealBins.data(),
            reducedAmpFactors.data(),
            fftSize);
    if (blockIndex == blockIndexToLog) {
        logSignal("100-input-fft-real-x-reduced-amp-factors-b.txt", fftSize, inputRealWithReducedAmpFactors.data());
    }

    // Multiply the imaginary components of the input fft by the reduced
    // combined amps.
    FFTArray inputImagWithReducedAmpFactors;
    FloatVectorOperations::multiply(
            inputImagWithReducedAmpFactors.data(),
            inputImagBins.data(),
            reducedAmpFactors.data(),
            fftSize);
    if (blockIndex == blockIndexToLog) {
        logSignal("150-input-fft-imag-x-reduced-amp-factors-b.txt", fftSize, inputImagWithReducedAmpFactors.data());
    }

    ComplexFFTArray ifftData;
    getIFFT(inputRealWithReducedAmpFactors, inputImagWithReducedAmpFactors, ifftData);
    //getIFFT(inputRealBins, inputImagBins, ifftData);

    applyHannWindow(ifftData.data(), ifftData.size());

    // Copy the results to the channel.
    for (int i = 0; i < fftSize; ++i) {
        outBlockSamples[i] = ifftData[i];
    }
}
