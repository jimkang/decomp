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
  int numberOfOuts, int blockIndex);

static void getReducedCombinedAmpFactors(ComplexFFTArray &inputFFTData, FFTArray &reducedAmpFactors);

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

    const int blockSize = 1 << numberOfOuts;

    const int maxBlocks = inputSamples.size() / blockSize;
    // Leave out the last partial block for now.

    auto inputStart = inputSamples.begin();

    IIRFilter inputHighPassFilter;
    double freq = 5;
    inputHighPassFilter.setCoefficients(IIRCoefficients::makeHighPass(sampleRate, freq));

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
  int numberOfOuts, int blockIndex) {

  const int blockSize = 1 << numberOfOuts;
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

  applyHannWindow(inputBlockSamples.data(), inputBlockSamples.size());

  // Run an FFT on the signal.
  std::array<float, fftSize * 2> inputFFTData;

  for (int i = 0; i < inputBlockSamples.size(); ++i) {
      inputFFTData[i] = inputBlockSamples[i];
  }

  const int fftOrder = numberOfOuts;
  const int fftSize = 1 << fftOrder;
  dsp::FFT fft(fftOrder);
  // TODO: Use the imaginary values, look up bandpass filter impl. for reference.
  fft.performRealOnlyForwardTransform(inputFFTData.data(), true);
  //getFFT(inputFFTData);

  auto maxLevel = juce::FloatVectorOperations::findMinAndMax (inputFFTData.data(), fftPowerOf2 / 2);

  for (int fftIndex = 0; fftIndex < fftSize/2; ++fftIndex) {
    const float freq = fftIndex * fftSize / 2;
    const float level = juce::jmap(
      inputFFTData[fftIndex],
      0.0f, juce::jmax(maxLevel.getEnd(), 1e-5f),
      0.0f, 1.0f
    );
    cout << "Block " << blockIndex << ": freq " << freq << ", level " << level << endl;
  }

    //applyHannWindow(ifftData.data(), ifftData.size());

    // Copy the results to the channel.
    //for (int i = 0; i < fftSize; ++i) {
    //   outBlockSamples[i] = ifftData[i];
    //}
  return outBlockSamplesCollection;
}
