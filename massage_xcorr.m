function c = massage_xcorr(c)
c = sort(c);
plot(c,'corr');
% set(gcf,'Position',[50 50 500 400]);
corr_matrix = get(c,'CORR');

plot(c,'lag');
% set(gcf,'Position',[50 50 500 400]);
lag_matrix = get(c,'LAG');
lag_matrix(1:5,1:5)

c = adjusttrig(c,'MIN',1);
plot(c,'sha');
% set(gcf,'Position',[50 50 700 400]);

c = linkage(c);
plot(c,'den');
close(gcf)
% set(gcf,'Position',[50 50 600 400]);

c = cluster(c,.8);
index = find(c,'CLUST',1);
c1 = subset(c,index);
plot(c1,'wig');
% set(gcf,'Position',[50 50 700 400]);

c1 = crop(c1,-4,9);
c1 = norm(c1);
c1 = stack(c1);
c1 = norm(c1);

c2 = norm(c1);
c2 = minus(c2);
plot(c2,'raw',.5)
% set(gcf,'Position',[50 50 700 400]);

c1 = xcorr(c1,[1 3]);       % align traces on the initial wavelet
c1 = adjusttrig(c1);        %     "     "
c1 = interferogram(c1);
plot(c1,'int',1,'corr');    % plot the correlation value interferogram
% set(gcf,'Position',[50 50 700 400]);
plot(c1,'int',1,'lag',.01);    % plot the lag value interferogram
% set(gcf,'Position',[50 50 700 400]);

c = crop(c,-3,10);
plot(c,'occurence',1,1:max(get(c,'CLUST')));
% set(gcf,'Position',[50 50 700 400]);

index = find(c,'CLUST',1);
plot(c,'overlay',1,index);
% set(gcf,'Position',[50 50 700 200]);
end