/*
  ==============================================================================

    This file contains the basic startup code for a JUCE application.

  ==============================================================================
*/

#include <JuceHeader.h>
#include "Decomp.h"
#include "Reconstruct.h"
#include "WriteAudio.h"
#include "RunStep.h"

static const char *const usageMsg = "Usage: decomp input.wav <freq. bucket count> <out dir>";

//==============================================================================
int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << usageMsg << std::endl;
        return 1;
    }

    juce::File inputFile(argv[1]);
    if (!inputFile.existsAsFile()) {
        std::cerr << "Input file does not exist." << std::endl;
        return 1;
    }

    const int bucketCount = atoi(argv[2]);

    std::string outDirPath = argv[3];
    juce::File outDir(outDirPath);
    if (!outDir.isDirectory()) {
        std::cerr << "Output directory does not exist." << std::endl;
        return 1;
    }

    juce::AudioFormatManager formatManager;
    formatManager.registerBasicFormats();

    juce::AudioFormatReader *inputReader = formatManager.createReaderFor(inputFile);
    if (inputReader == nullptr) {
        std::cerr << "Could not read input file." << std::endl;
        return 1;
    }

    juce::AudioBuffer<float> inputBuffer;
    inputBuffer.setSize(inputReader->numChannels, inputReader->lengthInSamples);
    inputReader->read(&inputBuffer, 0, inputReader->lengthInSamples, 0, true, true);
    std::cout << "Read " << inputReader->lengthInSamples << " samples from input." << std::endl;

    const double sampleRate = inputReader->sampleRate;
    const int bitsPerSample = inputReader->bitsPerSample;

    delete inputReader;

    const int outLen = inputReader->lengthInSamples;
    const int channelCount = inputReader->numChannels;

    vector<juce::AudioBuffer<float>> outBuffers;
    for (int i = 0; i < bucketCount; ++i) {
      juce::AudioBuffer<float> outBuffer;
      outBuffer.setSize(channelCount, outLen);
      outBuffers.push_back(outBuffer);
    }

    //reconstruct(inputBuffer, outBuffer);
    decomp(inputBuffer, sampleRate, bucketCount, outBuffers);

    int bufferNumber = 0;
    for (auto it = outBuffers.begin(); it != outBuffers.end(); ++it) {
      std::string outFilePath(outDirPath);
      outFilePath += "/";
      outFilePath += std::to_string(bufferNumber);
      outFilePath += ".wav";
      writeBufferToFile(*it, sampleRate, bitsPerSample, outFilePath.c_str());
      bufferNumber += 1;
    }

    return 0;
}
