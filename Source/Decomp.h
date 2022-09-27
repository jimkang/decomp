#include <JuceHeader.h>
#include <math.h>
#include "Utils.h"
#include "Consts.h"
#include "LogSignal.h"
#include "Sqrt.h"
#include "steps/HannWindowStep.h"

using namespace juce;
using namespace std;

static vector<vector<float>> decompChannel(
  vector<float> &inputSamples, double sampleRate, int numberOfOuts);

static vector<vector<float>> decompBlock(
  vector<float> &inputBlockSamples, IIRFilter &inputHighPassFilter,
  IIRFilter &infoHighPassFilter, int numberOfOuts, int blockIndex);

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
    vector<vector<float>> outSamplesCollection = decompChannel(
      inputSamples, sampleRate, outBuffers.size()
    );
    cout << "outSamplesCollection length: " << outSamplesCollection.size() << endl;

    for (int outNumber = 0; outNumber < outSamplesCollection.size(); ++outNumber) {
      cout << "outNumber: " << outNumber << endl;
      float *outWritePtr = outBuffers[outNumber].getWritePointer(ch);
      vector<float> outputSamples = outSamplesCollection[outNumber];
      cout << "outputSamples length at end: " << outputSamples.size() << endl;
      cout << "outputSamples value at end: " << outputSamples[1000] << endl;
      for (int i = 0; i < outputSamples.size(); ++i) {
        outWritePtr[i] = outputSamples[i];
      }
    }
  }
}

static vector<vector<float>> decompChannel(vector<float> &inputSamples,
  double sampleRate, int numberOfOuts) {

    const int maxBlocks = inputSamples.size() / blockSize;
    // Leave out the last partial block for now.

    auto inputStart = inputSamples.begin();

    IIRFilter inputHighPassFilter;
    IIRFilter infoHighPassFilter;
    double freq = 5;
    inputHighPassFilter.setCoefficients(IIRCoefficients::makeHighPass(sampleRate, freq));
    infoHighPassFilter.setCoefficients(IIRCoefficients::makeHighPass(sampleRate, freq));

    vector<vector<float>> outChannelSamplesCollection;
    for (int i = 0; i < numberOfOuts; ++i) {
      vector<float> outChannelSamples(inputSamples.size());
      outChannelSamplesCollection.push_back(outChannelSamples);
    }

    for (int blockIndex = 0; blockIndex < maxBlocks; ++blockIndex) {
      cout << "analyzing block  " << blockIndex << endl;
      const int sampleIndex = blockIndex * blockSize;

      vector<float> inputBlockSamples(blockSize);
      vector<vector<float>> outBlockSamplesCollection = decompBlock(
        inputBlockSamples,
        inputHighPassFilter,
        infoHighPassFilter,
        numberOfOuts,
        blockIndex
      );
      cout << "outBlockSamplesCollection length:" << outBlockSamplesCollection.size() << endl;
      cout << "outBlockSamplesCollection sample:" << outBlockSamplesCollection[3][1000] << endl;

      // Copy block results to the larger sample vector.
      // TODO: Cut down on the copying.
      for (int outNumber = 0; outNumber < outBlockSamplesCollection.size(); ++outNumber) {
        auto outBlockSamples = outBlockSamplesCollection[outNumber];
        auto outChannelSamples = outChannelSamplesCollection[outNumber];

        cout << "outBlockSamples sample:" << outBlockSamples[1000] << endl;
        //std::copy(outBlockSamples.begin(), outBlockSamples.end(), std::back_inserter(outChannelSamples));
        for (int i = 0; i < outBlockSamples.size(); ++i) {
          const int destIndex = blockSize * blockIndex + i;
          if (destIndex >= outChannelSamples.size()) {
            break;
          }
          outChannelSamples[destIndex] = outBlockSamples[blockIndex];
        }

        cout << "outBlockSamples length:" << outBlockSamples.size() << endl;
        cout << "outChannelSamples length:" << outChannelSamples.size() << endl;
        // TODO: This is necessary, so it's a sign the overall approach is wrong.
        outChannelSamplesCollection[outNumber] = outChannelSamples;
      }
    }
    cout << "End of decompChannel outChannelSamplesCollection length: " << outChannelSamplesCollection.size() << endl;
    cout << "End of decompChannel sample: " << outChannelSamplesCollection[3][1000] << endl;
    return outChannelSamplesCollection;
}

static vector<vector<float>> decompBlock(
  vector<float> &inputBlockSamples, IIRFilter &inputHighPassFilter,
  IIRFilter &infoHighPassFilter, int numberOfOuts, int blockIndex) {

  vector<vector<float>> outBlockSamplesCollection;

  auto inputSample = inputBlockSamples.begin();
  for (int outNumber = 0; outNumber < numberOfOuts; ++outNumber) {
    vector<float> outBlockSamples;
    for (int i = 0; i < blockSize; ++i) {
      outBlockSamples.push_back(1.0);//inputSample[i];
    }
    cout << "outBlockSamples.size()" << outBlockSamples.size() << endl;
    outBlockSamplesCollection.push_back(outBlockSamples);
    cout << "outBlockSamplesCollection.back().size()" << outBlockSamplesCollection.back().size() << endl;
    cout << "outBlockSamplesCollection.size()" << outBlockSamplesCollection.size() << endl;
  }

  return outBlockSamplesCollection;

   // Drop high-pass for now.
    //inputHighPassFilter.processSamples(inputBlockSamples.data(), inputBlockSamples.size());
    //infoHighPassFilter.processSamples(infoBlockSamples.data(), infoBlockSamples.size());

    applyHannWindow(inputBlockSamples.data(), inputBlockSamples.size());

    // Run a real-only FFT on both signals.
    ComplexFFTArray inputFFTData;

    for (int i = 0; i < inputBlockSamples.size(); ++i) {
        inputFFTData[i] = inputBlockSamples[i];
    }

    getFFT(inputFFTData);

    //if (blockIndex == blockIndexToLog) {
    //logSignal("040-inputFFT-b.txt", fftSize, inputFFTData.data());
    //}
/*
    FFTArray infoRealBins;
    getReal(infoFFTData, infoRealBins);

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
*/
}
