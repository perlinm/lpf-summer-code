function output = signal_filter(input,low_pass_filter,band_pass_filter)
    % act low pass filter
    output = fftfilt(input,low_pass_filter);
    % downsample to frequency of band pass filter
    output.downsample(plist(...
        'factor',input.fs/band_pass_filter.hist.plistUsed.find('fs')));
    % act band pass filter
    output.fftfilt(band_pass_filter);
    % remove transient data
    output.select(...
        band_pass_filter.hist.plistUsed.find('order'):len(output));
end