function wf = getPkFreq(wf,plotz)
%% Finds the maximum frequency for each element of the waveform wf
  % syntax: wf = getPkFreq(wf,false) or wf = getPkFreq(wf,true)
  %   Inputs: wf = waveform object
  %           plotz = true for a plot, false for no plot
  %   Output: wf = waveform object with added fields:
  %        FFT_FREQ = array of frequencies from the fft of the waveform
  %        data
  %        FFT_AMP = strength of the signal at each frequency in FFT_FREQ
  %        FFT_DOM = dominant frequency of the waveform data
    for count = 1:numel(wf)

        Fn = get(wf(count),'NYQ'); %from https://code.google.com/p/gismotools/source/browse/trunk/GISMO/contributed/fft_tools/%2Bwf_fft/compute.m?r=321
        stn = get(wf(count),'station');
        x = get(wf(count),'data');


        NFFT=2.^(ceil(log(length(x))/log(2)));  % Next highest power of 2
        FFTX=fft(x,NFFT);                       % Take fft, padding with zeros.
        NumUniquePts = ceil((NFFT+1)/2);
        FFTX=FFTX(1:NumUniquePts);              % throw out neg frequencies
        MX=abs(FFTX);                           % Take magnitude of X
        MX=MX*2;                                % Multiply by 2 
        MX=MX/length(x);                        
        %PX=phase(FFTX);                           % Take magnitude of X
        f=(0:NumUniquePts-1)*2/NFFT;            
        f=f*Fn;
        [~,I] = max(MX);

        if plotz
            figure
            plot(f,MX)
            title(sprintf('%s',stn)); hold on
            line([f(I), f(I)], [min(MX), max(MX)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
        end
        wf(count) = addfield(wf(count),'FFT_FREQ',f');
        wf(count) = addfield(wf(count),'FFT_AMP',MX);
        wf(count) = addfield(wf(count),'FFT_DOM',f(I));
        clearvars -except count wf plotz
    end
end