function pk_vals = det_pks(signals,frequencies,min_freq,type)  
% min_freq in Hz

% signal for periodogram needs no manipulation

%% Determine peaks for stft
    if strcmp(type, 'STFT')
        signals = mean(signals,2);
        
        
    end
    
    [PKS,LOCS] = findpeaks(signals);
    freq_inds = find(frequencies(LOCS)>=min_freq);
    PK_freq = frequencies(LOCS(freq_inds));
    best_fit_pks = find(PKS(freq_inds)>= mean(PKS(freq_inds))+range(PKS(freq_inds))/6);
    pk_vals = PK_freq(best_fit_pks);