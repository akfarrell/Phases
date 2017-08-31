%function short_time_freq_comp()

% Try different techniques to get the best short time window frequency information

    %This includes :
    % 1.) periodogram
    % 2.) spectrogram
    % 3.) short time fourier transform
    % 4.) multitaper
    
    %this works by creating two signals:
    % 1.) sine wave of two frequencies with random noise
    % 2.) same sine wave of two frequencies but with trend
    close all; clc;
    addpath('/home/a/akfarrell/Uturuncu/Phase/stft')
    indz = 5:8;
    sig_lengths = 2.^indz;
    analy = {'STFT','Periodogram'};
    
    %% Generate signals
    for count1 = 1:numel(sig_lengths)
        sig_len = sig_lengths(count1);
        x = 0:sig_len-1; %128(127) okay, 256(255) is best
        f1 = 12;
        f2 = 17;
        fs = 100;
        y = sin(2*pi*f1/fs*x)+sin(2*pi*f2/fs*x);
        sig1 = y+0.15*randn(1,length(y)); %adding 0.15 amplitude units of randnoise


        sig2 = sig1+x/(numel(x)/5);


        %% Fourier transform inputs
        leng = numel(x);
        wlen = leng/2; %x/2
        nfft = leng; %x
        hop = leng/4; %x/4

        %% Compare different methods
        [pdgrm,f_pdgrm] = periodogram(sig1,[],nfft,fs);
        [pdgrm_sig2,f_pdgrm_sig2] = periodogram(sig2,[],nfft,fs);
        %spcgrm = ;
        [stft_res,f_stft,t_stft] = stft(sig1, wlen, hop, nfft, fs);
        [stft_res_sig2,f_stft_sig2,t_stft_sig2] = stft(sig2, wlen, hop, nfft, fs);
        %mlttpr = ;

        %% Plot the signals and methods
        for count = 1:numel(analy)
            p =figure();
            set(p, 'Position', [1000 1000 1200 1200])
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)    
            subplot(2,2,1)
                plot(x*1/fs,sig1);xlim([min(x*1/fs) max(x*1/fs)])
                xlabel('Time (s)')
                ylabel('Amplitude')
                title(sprintf('Signal 1, %d to %d Hz & random noise %d samples',f1,f2,sig_len))
            subplot(2,2,2)
                plot(x*1/fs,sig2);xlim([min(x*1/fs) max(x*1/fs)])
                xlabel('Time (s)')
                ylabel('Amplitude')
                title(sprintf('Signal 2, same as Signal 1 but with trend %d samples',sig_len))
            subplot(2,2,3)
            if strcmp(analy{count},'STFT')
                    stft_res_sig1_db = stft_plotting(wlen, f_stft,nfft,stft_res,t_stft);
                subplot(2,2,4)
                    stft_res_sig2_db = stft_plotting(wlen, f_stft_sig2,nfft,stft_res_sig2,t_stft_sig2);
            elseif strcmp(analy{count},'Periodogram')
                    plot(f_pdgrm,pdgrm)
                    xlabel('Frequency (Hz)')
                    ylabel('Power Spectral Density')
                    title(sprintf('%s Signal 1',analy{count}))
                subplot(2,2,4)
                    plot(f_pdgrm_sig2,pdgrm_sig2)
                    xlabel('Frequency (Hz)')
                    ylabel('Power Spectral Density')
                    title(sprintf('%s Signal 2',analy{count}))
            end
            %% Saving file
                hold off
                directory = '/home/a/akfarrell/Uturuncu/Phase/comparison';
                filename = sprintf('compare_%s_%dsamples',analy{count},leng);
                filename_wPath = fullfile(directory,filename);
                hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
        end

        %% Determine dominant frequency peaks
        pk_vals.(sprintf('length_%d',leng)).sig1.stft = det_pks(stft_res_sig1_db,f_stft,2,analy{1});
        pk_vals.(sprintf('length_%d',leng)).sig2.stft = det_pks(stft_res_sig2_db,f_stft_sig2,2,analy{1});
        pk_vals.(sprintf('length_%d',leng)).sig1.pdgrm = det_pks(pdgrm,f_pdgrm,2,analy{2})';
        pk_vals.(sprintf('length_%d',leng)).sig2.pdgrm = det_pks(pdgrm_sig2,f_pdgrm_sig2,2,analy{2})';
        clearvars -except pk_vals count1 analy sig_lengths
    end
