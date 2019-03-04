close all
load('wf_objs/wf_2166.mat');
%s = spectrogram(get(w_clean(1),'data'), floor(numel(get(w_clean(1),'data'))/500), floor(numel(get(w_clean(1),'data'))/500)-1, 128, 100);
figure()
subplot(1,2,1)
plot(get(w_clean(1),'data'),1:numel(get(w_clean(1),'data')))
ylim([1 numel(get(w_clean(1),'data'))])
hold on
subplot(1,2,2)
spectrogram(get(w_clean(1),'data'),31,30,128,100);
