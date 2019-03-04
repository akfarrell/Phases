%function short_time_freq_comp_iterations()

% Try different techniques to get the best short time window frequency information

    %This includes :
    % 1.) periodogram
    % 2.) spectrogram
    % 3.) short time fourier transform
    % 4.) multitaper
    
    %this works by creating two signals:
    % 1.) sine wave of two frequencies with random noise
    % 2.) same sine wave of two frequencies but with trend
    close all; clc; clear all;
    addpath('/home/a/akfarrell/Uturuncu/Phase/stft')
    indz = 5:8;
    sig_lengths = 2.^indz;
    analy = {'STFT','Periodogram','Spectrogram','Multitaper'};
    freq_cutoff = 2;
    num_f = 2; %CHANGE!!!
    f1 = 12;
    %f2 = 13;
    f2 = 17;
    %f3 = 16;
    %f4  =16;
    %f5 = 17;
    stft_sig1 = zeros(100,num_f+1,4);
    pdgrm_sig1 = zeros(100,num_f+1,4);
    spcgrm_sig1 = zeros(100,num_f+1,4);
    mlttpr_sig1 = zeros(100,num_f+1,4);
    
    stft_sig2 = zeros(100,num_f+1,4);
    pdgrm_sig2 = zeros(100,num_f+1,4);
    spcgrm_sig2 = zeros(100,num_f+1,4);
    mlttpr_sig2 = zeros(100,num_f+1,4);
    
    for count = 1:100
        %% Generate signals
        for count1 = 1:numel(sig_lengths)
            sig_len = sig_lengths(count1);
            x = 0:sig_len-1; %128(127) okay, 256(255) is best
            fs = 100;
            y = sin(2*pi*f1/fs*x)+0.5*sin(2*pi*f2/fs*x);%+sin(2*pi*f3/fs*x);%+sin(2*pi*f4/fs*x);%+sin(2*pi*f5/fs*x);
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
            [spcgrm, f_spcgrm,t_spcgrm] = spectrogram(sig1,wlen,wlen/2,nfft,fs);
            [spcgrm_sig2, f_spcgrm_sig2,t_spcgrm_sig2] = spectrogram(sig2,wlen,wlen/2,nfft,fs);
            [stft_res,f_stft,t_stft] = stft(sig1, wlen, hop, nfft, fs);
            [stft_res_sig2,f_stft_sig2,t_stft_sig2] = stft(sig2, wlen, hop, nfft, fs);
            [mlttpr, f_mlttpr] = pmtm(sig1,2,nfft,fs,'unity');
            [mlttpr_sig2, f_mlttpr_sig2] = pmtm(sig2,2,nfft,fs,'unity');
            stft_res_sig1_db = stft_plotting(wlen, f_stft,nfft,stft_res,t_stft,'STFT');
            stft_res_sig2_db = stft_plotting(wlen, f_stft_sig2,nfft,stft_res_sig2,t_stft_sig2,'STFT');
    %         %% Try using dpss ------ May not work!!!! ------

    % Doesn't improve values for multitaper; multitaper still worst method

    %         [e1,v1] = dpss(length(sig1),2,10);
    %         idx1 = find(v1>0.99,1,'last');
    %         [pxx1,f3] = pmtm(sig1,e1(:,1:idx1),v1(1:idx1),length(sig1),fs,'DropLastTaper',false);
    %         det_pks(pxx1,f3,2,'Multitaper')
    %         
    %         [e2,v2] = dpss(length(sig2),2,10);
    %         idx2 = find(v2>0.99,1,'last');
    %         [pxx2,f4] = pmtm(sig2,e2(:,1:idx2),v2(1:idx2),length(sig2),fs,'DropLastTaper',false);
    %         det_pks(pxx2,f4,2,'Multitaper')

%             end
%% Save values to arrays
if count1 == 1
    k = det_pks(stft_res_sig1_db,f_stft,freq_cutoff,analy{1});
    stft_sig1(count,1:numel(k)) = k;
    clear k; k=det_pks(pdgrm,f_pdgrm,freq_cutoff,analy{2})';
    pdgrm_sig1(count,1:numel(k)) = k;
    clear k; k=det_pks(abs(spcgrm),f_spcgrm,freq_cutoff,analy{3})';
    spcgrm_sig1(count,1:numel(k)) = k;
    clear k; k=det_pks(mlttpr,f_mlttpr,freq_cutoff,analy{4})';
    mlttpr_sig1(count,1:numel(k)) = k;
    
    clear k; k=det_pks(stft_res_sig2_db,f_stft_sig2,freq_cutoff,analy{1});
    stft_sig2(count,1:numel(k)) = k;
    clear k; k=det_pks(pdgrm_sig2,f_pdgrm_sig2,freq_cutoff,analy{2})';
    pdgrm_sig2(count,1:numel(k)) = k;
    clear k; k=det_pks(abs(spcgrm_sig2),f_spcgrm_sig2,freq_cutoff,analy{3})';
    spcgrm_sig2(count,1:numel(k)) = k;
    clear k; k=det_pks(mlttpr_sig2,f_mlttpr_sig2,freq_cutoff,analy{4})'; 
    mlttpr_sig2(count,1:numel(k)) =   k; 
elseif count1 == 2
    clear k; k=det_pks(stft_res_sig1_db,f_stft,freq_cutoff,analy{1});
    stft_sig1(count,1:numel(k),2) = k;
    clear k; k=det_pks(pdgrm,f_pdgrm,freq_cutoff,analy{2})';
    pdgrm_sig1(count,1:numel(k),2) = k;
    clear k; k=det_pks(abs(spcgrm),f_spcgrm,freq_cutoff,analy{3})';
    spcgrm_sig1(count,1:numel(k),2) = k;
    clear k; k=det_pks(mlttpr,f_mlttpr,freq_cutoff,analy{4})';
    mlttpr_sig1(count,1:numel(k),2) = k;
    
    clear k; k=det_pks(stft_res_sig2_db,f_stft_sig2,freq_cutoff,analy{1});
    stft_sig2(count,1:numel(k),2) = k;
    clear k; k=det_pks(pdgrm_sig2,f_pdgrm_sig2,freq_cutoff,analy{2})';
    pdgrm_sig2(count,1:numel(k),2) = k;
    clear k; k=det_pks(abs(spcgrm_sig2),f_spcgrm_sig2,freq_cutoff,analy{3})';
    spcgrm_sig2(count,1:numel(k),2) = k;
    clear k; k=det_pks(mlttpr_sig2,f_mlttpr_sig2,freq_cutoff,analy{4})';
    mlttpr_sig2(count,1:numel(k),2) =  k;
elseif count1 == 3
    clear k; k=det_pks(stft_res_sig1_db,f_stft,freq_cutoff,analy{1});
    stft_sig1(count,1:numel(k),3) = k;
    clear k; k=det_pks(pdgrm,f_pdgrm,freq_cutoff,analy{2})';
    pdgrm_sig1(count,1:numel(k),3) = k;
    clear k; k=det_pks(abs(spcgrm),f_spcgrm,freq_cutoff,analy{3})';
    spcgrm_sig1(count,1:numel(k),3) = k;
    clear k; k=det_pks(mlttpr,f_mlttpr,freq_cutoff,analy{4})';
    mlttpr_sig1(count,1:numel(k),3) = k;
    
    clear k; k=det_pks(stft_res_sig2_db,f_stft_sig2,freq_cutoff,analy{1});
    stft_sig2(count,1:numel(k),3) = k;
    clear k; k=det_pks(pdgrm_sig2,f_pdgrm_sig2,freq_cutoff,analy{2})';
    pdgrm_sig2(count,1:numel(k),3) = k;
    clear k; k=det_pks(abs(spcgrm_sig2),f_spcgrm_sig2,freq_cutoff,analy{3})';
    spcgrm_sig2(count,1:numel(k),3) = k;
    clear k; k=det_pks(mlttpr_sig2,f_mlttpr_sig2,freq_cutoff,analy{4})';  
    mlttpr_sig2(count,1:numel(k),3) = k;
elseif count1 == 4
    clear k; k=det_pks(stft_res_sig1_db,f_stft,freq_cutoff,analy{1});
    stft_sig1(count,1:numel(k),4) = k;
    clear k; k=det_pks(pdgrm,f_pdgrm,freq_cutoff,analy{2})';
    pdgrm_sig1(count,1:numel(k),4) = k;
    clear k; k=det_pks(abs(spcgrm),f_spcgrm,freq_cutoff,analy{3})';
    spcgrm_sig1(count,1:numel(k),4) = k;
    clear k; k=det_pks(mlttpr,f_mlttpr,freq_cutoff,analy{4})';
    mlttpr_sig1(count,1:numel(k),4) = k;
    
    clear k; k=det_pks(stft_res_sig2_db,f_stft_sig2,freq_cutoff,analy{1});
    stft_sig2(count,1:numel(k),4) = k;
    clear k; k=det_pks(pdgrm_sig2,f_pdgrm_sig2,freq_cutoff,analy{2})';
    pdgrm_sig2(count,1:numel(k),4) = k;
    clear k; k=det_pks(abs(spcgrm_sig2),f_spcgrm_sig2,freq_cutoff,analy{3})';
    spcgrm_sig2(count,1:numel(k),4) = k;
    clear k; k=det_pks(mlttpr_sig2,f_mlttpr_sig2,freq_cutoff,analy{4})'; 
    mlttpr_sig2(count,1:numel(k),4) =  k;
end


            clearvars -except pk_vals count1 analy sig_lengths freq_cutoff f1 f2 f3 f4 f5 num_f mlttpr_s* spcgrm_s* stft_s* pdgrm_s* count
        end

    end
%[a,b] = hist(mlttpr_sig1(:,1,1),unique(mlttpr_sig1(:,1,1)))
%     for count2 = 1:2 %signal 1 and signal 2
%         for count3 = 1:numel(num_f)
%             for count4 = 1:numel(analy)
%                 for count5 = 1:numel(sig_lengths)
%                     
%                 end
%             end
%         end
%     end