function [pk_vals,width,which_pk] = det_pks(signals,frequencies,min_freq,type,param)  
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
    
    [PKS,LOCS,W] = findpeaks(signals,frequencies)
    freq_inds = find(LOCS>=min_freq);
    ind = find(PKS(freq_inds)>= mean(PKS(freq_inds))+range(PKS(freq_inds))/6)
    if numel(ind) > 1
        %PKS(freq_inds(ind));
        PKS(freq_inds(ind))
        ind = find(PKS(freq_inds) == max(PKS(freq_inds(ind))))
        
    end
    if exist('param','var')
        if ind >1
            (PKS(freq_inds(ind-1))/PKS(freq_inds(ind))<(3/4))
            PKS(ind-1)/PKS(ind)
            if ~(PKS(freq_inds(ind-1))/PKS(freq_inds(ind))<(3/4))
                ind = ind-1;
                which_pk = '-';
            else
                disp('Lower peak too low in amplitude to consider!!!!!')
                which_pk = '+';
            end
            pk_vals = LOCS(freq_inds(ind));
            width = W(ind);
        else
            pk_vals = LOCS(freq_inds(ind));
            width = W(ind);
            disp('No lower peak - had to go with higher peak!!!')
            which_pk = '+';
        end
    end
    
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