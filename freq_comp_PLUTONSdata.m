%function freq_comp_PLUTONSdata()
  
    %this works by creating two signals:
    % 1.) sine wave of two frequencies with random noise
    % 2.) same sine wave of two frequencies but with trend
    %close all; clc; clear all;
    addpath('/home/a/akfarrell/Uturuncu/Phase/stft')
    %indz = 5:8;
    %sig_lengths = 2.^indz;
    analy = {'STFT','Periodogram','Spectrogram','Multitaper'};
    start_value = 547; %starting value in Haystack_data of the phase
    
    indexes = start_value:start_value+31;
    freq_cutoff = 2;
    num_f = 2; %CHANGE!!!
    fs=100;
    
    stft_sig = zeros(1,num_f+1,3);
    pdgrm_sig = zeros(1,num_f+1,3);
    spcgrm_sig = zeros(1,num_f+1,3);
    mlttpr_sig = zeros(1,num_f+1,3);
    
    sigz(1:32,1,1) = Haystack_data.HHZ(indexes);
    sigz(1:32,1,2) = Haystack_data.HHR(indexes);
    sigz(1:32,1,3) = Haystack_data.HHT(indexes);
    
    for count1 = 1:3
        sig = sigz(:,1,count1);
      
        %% Fourier transform inputs
        leng = numel(sig); %%%%%%%%CHANGE TO NUMBER OF POINTS IN SIGNAL!!!!!!
        wlen = leng/2; %x/2
        nfft = leng; %x
        hop = leng/4; %x/4

        %% Compare different methods
        [pdgrm,f_pdgrm] = periodogram(sig,[],nfft,fs);
        [Peak,PeakIdx]=findpeaks(pdgrm,f_pdgrm);
        figure()
        plot(f_pdgrm,pdgrm)
        text(PeakIdx,Peak+5, num2cell(PeakIdx,2))
        [spcgrm, f_spcgrm,t_spcgrm] = spectrogram(sig,wlen,wlen/2,nfft,fs);
        [stft_res,f_stft,t_stft] = stft(sig, wlen, hop, nfft, fs);
        [mlttpr, f_mlttpr] = pmtm(sig,2,nfft,fs,'unity');
        stft_res_sig_db = stft_plotting(wlen, f_stft,nfft,stft_res,t_stft,'STFT');

        %% Save values to arrays
        k = det_pks(stft_res_sig_db,f_stft,freq_cutoff,analy{1});
        stft_sig(1,1:numel(k),count1) = k;
        clear k; k=det_pks(pdgrm,f_pdgrm,freq_cutoff,analy{2})';
        pdgrm_sig(1,1:numel(k),count1) = k;
        clear k; k=det_pks(abs(spcgrm),f_spcgrm,freq_cutoff,analy{3})';
        spcgrm_sig(1,1:numel(k),count1) = k;
        clear k; k=det_pks(mlttpr,f_mlttpr,freq_cutoff,analy{4})';
        mlttpr_sig(1,1:numel(k),count1) = k;

        %clearvars -except analy freq_cutoff num_f mlttpr_sig spcgrm_sig stft_sig pdgrm_sig count sigz count1 sig fs Haystack_data

    end