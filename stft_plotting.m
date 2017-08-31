function stft_res_db = stft_plotting(wlen, f_stft,nfft,stft_res,t_stft)
%% STFT plotting as part of a pre-created figure or subplot
    K = sum(hamming(wlen, 'periodic'))/wlen;
    stft_res = abs(stft_res)/wlen/K;
    if rem(nfft, 2)                     % odd nfft excludes Nyquist point
        stft_res(2:end, :) = stft_res(2:end, :).*2;
    else                                % even nfft includes Nyquist point
        stft_res(2:end-1, :) = stft_res(2:end-1, :).*2;
    end
    stft_res = 20*log10(stft_res + 1e-6); %convert to dB
    surf(t_stft, f_stft, stft_res)
    shading interp
    axis tight
    box on
    view(0, 90)
    xlabel('Time, s')
    ylabel('Frequency, Hz')
    title('Amplitude spectrogram of the signal from STFT')

    handl = colorbar;
    ylabel(handl, 'Magnitude, dB')
    stft_res_db = stft_res;
 