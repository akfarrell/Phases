figure()
subplot(4,1,1)
plot(lags,abs(c.HHE))
xlim([lags(1) lags(end)])
hold on
line([lags(1) lags(end)],[cutoff_val cutoff_val],'Color','r')
title('xcorr HHE')

subplot(4,1,2)
plot(1:numel(needle.HHE),needle.HHE)
title('Needle HHE')
xlim([1 numel(needle.HHE)])

subplot(4,1,3)
plot(1:numel(Haystack_data.HHE),Haystack_data.HHE)
title('Haystack HHE')
xlim([1 numel(Haystack_data.HHE)])
hold on 
plot(i.HHE:i.HHE+numel(needle.HHE)-1, needle.HHE,'r')

subplot(4,1,4)
plot(1:numel(needle.HHE),needle.HHE-Haystack_data.HHE(i.HHE:i.HHE+numel(needle.HHE)-1))
xlim([1 numel(needle.HHE)])
title('Diff Needle and Haystack Maximum Correlation')

figure()
subplot(4,1,1)
plot(lags,abs(c.HHN))
xlim([lags(1) lags(end)])
hold on
line([lags(1) lags(end)],[cutoff_val cutoff_val],'Color','r')
title('xcorr HHN')

subplot(4,1,2)
plot(1:numel(needle.HHN),needle.HHN)
title('Needle HHN')
xlim([1 numel(needle.HHN)])

subplot(4,1,3)
plot(1:numel(Haystack_data.HHN),Haystack_data.HHN)
title('Haystack HHN')
xlim([1 numel(Haystack_data.HHN)])
hold on 
plot(i.HHN:i.HHN+numel(needle.HHN)-1, needle.HHN,'r')

subplot(4,1,4)
plot(1:numel(needle.HHN),needle.HHN-Haystack_data.HHN(i.HHN:i.HHN+numel(needle.HHN)-1))
xlim([1 numel(needle.HHN)])
title('Diff Needle and Haystack Maximum Correlation')

figure()
subplot(4,1,1)
plot(lags,abs(c.HHZ))
xlim([lags(1) lags(end)])
hold on
line([lags(1) lags(end)],[cutoff_val cutoff_val],'Color','r')
title('xcorr HHZ')

subplot(4,1,2)
plot(1:numel(needle.HHZ),needle.HHZ)
title('Needle HHZ')
xlim([1 numel(needle.HHZ)])

subplot(4,1,3)
plot(1:numel(Haystack_data.HHZ),Haystack_data.HHZ)
title('Haystack HHZ')
xlim([1 numel(Haystack_data.HHZ)])
hold on 
plot(i.HHZ:i.HHZ+numel(needle.HHZ)-1, needle.HHZ,'r')

subplot(4,1,4)
plot(1:numel(needle.HHZ),needle.HHZ-Haystack_data.HHZ(i.HHZ:i.HHZ+numel(needle.HHZ)-1))
xlim([1 numel(needle.HHZ)])
title('Diff Needle and Haystack Maximum Correlation')
% end