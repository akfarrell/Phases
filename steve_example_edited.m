%% define datasources, channels etc
minutes_in_day = 60 * 24;
dbpath = '/raid/data/sakurajima/db';
ds = datasource('antelope', dbpath);
ctag(1) = ChannelTag('JP.SAKB.-.HHZ');
ctag(2) = ChannelTag('JP.SAKB.-.BD1');

%% plot multiple Sakurajima events from Jun 7th 2015
eventday = '7-Jun-2015';
eventtime = {'02:16:30'; '03:08:00'; '05:44:50'; '03:30:30'; '02:02:40'; '03:38:00'; '09:36:10'; '04:57:40'};
w = [];
for eventnum = 1:numel(eventtime)
    snum = datenum(sprintf('%s %s', eventday, eventtime{eventnum}));
    enum = snum + 1/minutes_in_day;
    w0 = waveform(ds, ctag, snum, enum);
    w = [w; w0];
end
plot_panels(w,true)

%% change the aspect ratio of the plot
ah = get(gcf,'Children'); % get all axes (panels) that belong to this figure
xpos = 0.4; xwidth=0.2; %CHANGE THESE TO CHANGE ASPECT RATIO
for c=2:length(ah) % loop over axes
    oldpos = get(ah(c), 'Position'); % get current position
    set(ah(c), 'Position', [xpos, oldpos(2), xwidth, oldpos(4)]); % change position
end