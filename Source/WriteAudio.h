#include <JuceHeader.h>

// Returns true it it was able to write.
bool writeBufferToFile(const juce::AudioBuffer<float> &outBuffer,
                       float sampleRate, float bitsPerSample, const char *filePath) {

    juce::File outFile(filePath);
    outFile.deleteFile();
    // Creating an outStream with the FileOutputStream on the stack creates
    // a problem because the AudioFormatWriter expects to be able to delete
    // it when it is destroyed. So, when we delete the writer below, it
    // will try to delete the outStream, which it cannot!
    // https://github.com/juce-framework/JUCE/blob/master/modules/juce_audio_formats/format/juce_AudioFormatWriter.cpp#L61
    //juce::FileOutputStream outStream(outFile, outLen * channelCount);
    std::unique_ptr<juce::FileOutputStream> outStream = outFile.createOutputStream();
    if (outStream->failedToOpen()) {
        std::cerr << "Could not open output file." << std::endl;
        return false;
    }

    juce::WavAudioFormat wavFormat;
    juce::AudioFormatWriter *writer = wavFormat.createWriterFor(
            outStream.get(), sampleRate, outBuffer.getNumChannels(), bitsPerSample, {}, 0
    );
    if (writer == nullptr) {
        std::cerr << "Could not write to output file." << std::endl;
        return false;
    }
    outStream.release();
    writer->writeFromAudioSampleBuffer(outBuffer, 0, outBuffer.getNumSamples());

    delete writer;

    return true;
}

// Returns true it it was able to write.
bool writeArrayToFile(const float *const *data, int channelCount, int sampleCount,
                      float sampleRate, float bitsPerSample, const char *filePath) {

    juce::File outFile(filePath);
    outFile.deleteFile();

    std::unique_ptr<juce::FileOutputStream> outStream = outFile.createOutputStream();
    if (outStream->failedToOpen()) {
        std::cerr << "Could not open output file." << std::endl;
        return false;
    }

    juce::WavAudioFormat wavFormat;
    juce::AudioFormatWriter *writer = wavFormat.createWriterFor(
            outStream.get(), sampleRate, channelCount, bitsPerSample, {}, 0
    );
    if (writer == nullptr) {
        std::cerr << "Could not write to output file." << std::endl;
        return false;
    }
    outStream.release();
    writer->writeFromFloatArrays(data, channelCount, sampleCount);

    delete writer;

    return true;
}
