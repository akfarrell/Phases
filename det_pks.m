function [pk_vals,width] = det_pks(signals,frequencies,min_freq,type)  
% min_freq in Hz

% signal for periodogram needs no manipulation
% signals = pdgrm1.(compnt{countz});
% frequencies = f_pdgrm1.(compnt{countz})
% min_freq = 3;
% type = 'pdgrm';


%% Determine peaks for stft
    if strcmp(type, 'STFT') || strcmp(type, 'Spectrogram')
        signals = mean(signals,2);
    end
    
    [PKS,LOCS,W] = findpeaks(signals,frequencies);
    freq_inds = find(LOCS>=min_freq);
    ind = find(PKS(freq_inds)>= mean(PKS(freq_inds))+range(PKS(freq_inds))/6)
    pk_vals = LOCS(freq_inds(ind));
    width = W(ind);
    
%     freq_inds = find(frequencies(LOCS)>=min_freq);
%     PK_freq = frequencies(LOCS(freq_inds));
%     best_fit_pks = find(PKS(freq_inds)>= mean(PKS(freq_inds))+range(PKS(freq_inds))/6);
%     pk_vals = PK_freq(best_fit_pks);
% ind = find(PKS(freq_inds)>= mean(PKS(freq_inds))+range(PKS(freq_inds))/6)
% pk_vals = PK_freq(ind);
% width = W(ind);
%%
%uncomment for figure!
figure()
findpeaks(signals,frequencies,'Annotate','extents')