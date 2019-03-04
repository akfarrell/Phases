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
    analy = {'STFT','Periodogram','Spectrogram','Multitaper'};
    freq_cutoff = 2;
    f1 = 12;
    f2 = 17;
    f3  =15;

    %% Generate signals
        for count1 = 1:numel(sig_lengths)
            sig_len = sig_lengths(count1);
            x = 0:sig_len-1; %128(127) okay, 256(255) is best
            fs = 100;
            y = sin(2*pi*f1/fs*x)+sin(2*pi*f2/fs*x)+sin(2*pi*f3/fs*x);
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

            %% Plot the signals and methods
             for count = 1:numel(analy)
%                 p =figure();
%                 set(p, 'Position', [1000 1000 1200 1200])
%                 set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)    
%                 subplot(2,2,1)
%                     plot(x*1/fs,sig1);xlim([min(x*1/fs) max(x*1/fs)])
%                     xlabel('Time (s)')
%                     ylabel('Amplitude')
%                     title(sprintf('Signal 1, %d to %d Hz & random noise %d samples',f1,f2,sig_len))
%                 subplot(2,2,2)
%                     plot(x*1/fs,sig2);xlim([min(x*1/fs) max(x*1/fs)])
%                     xlabel('Time (s)')
%                     ylabel('Amplitude')
%                     title(sprintf('Signal 2, same as Signal 1 but with trend %d samples',sig_len))
%                 subplot(2,2,3)
%                 if strcmp(analy{count},'STFT')
                        stft_res_sig1_db = stft_plotting(wlen, f_stft,nfft,stft_res,t_stft,analy{count});
%                     subplot(2,2,4)
                        stft_res_sig2_db = stft_plotting(wlen, f_stft_sig2,nfft,stft_res_sig2,t_stft_sig2,analy{count});
%                 elseif strcmp(analy{count},'Periodogram')
%                         plot(f_pdgrm,pdgrm)
%                         xlabel('Frequency (Hz)')
%                         ylabel('Power Spectral Density')
%                         title(sprintf('%s Signal 1',analy{count}))
%                     subplot(2,2,4)
%                         plot(f_pdgrm_sig2,pdgrm_sig2)
%                         xlabel('Frequency (Hz)')
%                         ylabel('Power Spectral Density')
%                         title(sprintf('%s Signal 2',analy{count}))
%                 elseif strcmp(analy{count},'Spectrogram')
%                         spec_res_sig1_db = stft_plotting(wlen,f_spcgrm,nfft,spcgrm,t_spcgrm,analy{count});
%                     subplot(2,2,4)
%                         spec_res_sig2_db = stft_plotting(wlen,f_spcgrm_sig2,nfft,spcgrm_sig2,t_spcgrm_sig2,analy{count});
%                 elseif strcmp(analy{count},'Multitaper')
%                         plot(f_mlttpr,mlttpr)
%                         xlabel('Frequency (Hz)')
%                         ylabel('Power Spectral Density')
%                         title(sprintf('%s Signal 1',analy{count}))
%                     subplot(2,2,4)
%                         plot(f_mlttpr_sig2,mlttpr_sig2)
%                         xlabel('Frequency (Hz)')
%                         ylabel('Power Spectral Density')
%                         title(sprintf('%s Signal 2',analy{count}))
                end
% 
%                 %% Saving file
%                     hold off
%                     directory = '/home/a/akfarrell/Uturuncu/Phase/comparison';
%                     filename = sprintf('compare_%s_%dsamples',analy{count},leng);
%                     filename_wPath = fullfile(directory,filename);
%                     hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

                %% Building another figure - all analysis methods for each sample number
%                 g = figure();
%                 set(g, 'Position', [1000 1000 1200 1200])
%                 set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)    
%                 subplot(2,2,1)
%                     plot(f_mlttpr,mlttpr)
%                     xlabel('Frequency (Hz)')
%                     ylabel('Power Spectral Density')
%                     title(sprintf('%s Signal 1',analy{4}))
%                 subplot(2,2,2)
%                     plot(f_pdgrm,pdgrm)
%                     xlabel('Frequency (Hz)')
%                     ylabel('Power Spectral Density')
%                     title(sprintf('%s Signal 1',analy{2}))
%                 subplot(2,2,3)
%                     stft_plotting(wlen,f_spcgrm,nfft,spcgrm,t_spcgrm,analy{3});
%                 subplot(2,2,4)    
%                     stft_plotting(wlen, f_stft,nfft,stft_res,t_stft,analy{1});
% 
%                 %% Saving file
%                 hold off
%                 directory = '/home/a/akfarrell/Uturuncu/Phase/comparison';
%                 filename = sprintf('compare_all_%dsamples_sig1',leng);
%                 filename_wPath = fullfile(directory,filename);
%                 hgexport(g, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%             end



%             %% Determine dominant frequency peaks
det_pks(stft_res_sig1_db,f_stft,freq_cutoff,analy{1})
%             pk_vals.(sprintf('length_%d',leng)).sig1.stft = det_pks(stft_res_sig1_db,f_stft,freq_cutoff,analy{1});
%             pk_vals.(sprintf('length_%d',leng)).sig2.stft = det_pks(stft_res_sig2_db,f_stft_sig2,freq_cutoff,analy{1});
%             pk_vals.(sprintf('length_%d',leng)).sig1.pdgrm = det_pks(pdgrm,f_pdgrm,freq_cutoff,analy{2})';
%             pk_vals.(sprintf('length_%d',leng)).sig2.pdgrm = det_pks(pdgrm_sig2,f_pdgrm_sig2,freq_cutoff,analy{2})';
%             pk_vals.(sprintf('length_%d',leng)).sig1.spcgrm = det_pks(abs(spcgrm),f_spcgrm,freq_cutoff,analy{3})';
%             pk_vals.(sprintf('length_%d',leng)).sig2.spcgrm = det_pks(abs(spcgrm_sig2),f_spcgrm_sig2,freq_cutoff,analy{3})';
%             pk_vals.(sprintf('length_%d',leng)).sig1.mlttpr = det_pks(mlttpr,f_mlttpr,freq_cutoff,analy{4})';
%             pk_vals.(sprintf('length_%d',leng)).sig2.mlttpr = det_pks(mlttpr_sig2,f_mlttpr_sig2,freq_cutoff,analy{4})';
            clearvars -except pk_vals count1 analy sig_lengths freq_cutoff f1 f2 f3 mlttpr_* spcgrm_* stft_* pdgrm_*
        end

        %% Use pk_vals to determine which method is consistently closer to the given frequencies
%         namez = fieldnames(pk_vals);
%         colorz = {'r','k'};
%         k = figure();
%         set(k, 'Position', [1000 1000 1200 1200])
%         set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
%         analy_type = {'stft','pdgrm','spcgrm','mlttpr'};
%         for count = 1:numel(namez)
%             plot_vals_sig1 = zeros(numel(analy_type),2);
%             plot_vals_sig2 = plot_vals_sig1;
%             for count2 = 1:numel(analy_type)
%                 plot_vals_sig1(count2,1:2) = pk_vals.(namez{count}).sig1.(analy_type{count2});
%                 plot_vals_sig2(count2,1:2) = pk_vals.(namez{count}).sig2.(analy_type{count2});
%             end
%             plot_num_s1 = [abs(plot_vals_sig1(:,1)-f1), abs(plot_vals_sig1(:,2)-f2)];
%             plot_num_s2 = [abs(plot_vals_sig2(:,1)-f1), abs(plot_vals_sig2(:,2)-f2)];
%             ax(1) = subplot(2,2,1);
%                 scatter(count,plot_num_s1(1,1),colorz{1},'f')
%                 hold on
%                 scatter(count,plot_num_s1(1,2),colorz{2},'f')
%                 if count == numel(namez)
%                     title('STFT')
%                 end
%             ax(2) = subplot(2,2,2);
%                 scatter(count,plot_num_s2(2,1),colorz{1},'f')
%                 hold on
%                 scatter(count,plot_num_s2(2,2),colorz{2},'f')
%                 if count == numel(namez)
%                     title('Periodogram')
%                 end
%             ax(3) = subplot(2,2,3);
%                 scatter(count,plot_num_s2(3,1),colorz{1},'f')
%                 hold on
%                 scatter(count,plot_num_s2(3,2),colorz{2},'f')
%                 if count == numel(namez)
%                     title('Spectrogram')
%                 end
%             ax(4) = subplot(2,2,4);
%                 scatter(count,plot_num_s2(4,1),colorz{1},'f')
%                 hold on
%                 scatter(count,plot_num_s2(4,2),colorz{2},'f')
%                 if count == numel(namez)
%                     title('Multitaper')
%                 end
%             clear count_array plot_num_s1 plot_num_s2 plot_vals_sig1 plot_vals_sig2
%         end            
%         linkaxes(ax,'xy')
