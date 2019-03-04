function phase_mulplt(w, alignWaveforms, max_vals, min_vals,fil,OrS,timez,oridStruct,style,station,subsetz)
%added last input for eq_noSwaves to tell difference between split-up
%waveform objects
%MULPLT Plot multiple waveform objects in a figure. is inspired by the 
%Seisan program of the same name
%   mulplt(w, alignWaveforms) 
%   where:
%       w = a vector of waveform objects
%       alignWaveforms is either true or false (default)
%   mulplt(w) will plot a record section, i.e. each waveform is plotted
%   against absolute time.
%   mulplt(w, true) will align the waveforms on their start times.

% Glenn Thompson 2014/11/05, generalized after a function I wrote in 2000
% to operate on Seisan files only
% Modified 4/2015 by Alexandra Farrell

    %w = waveform_nonempty(w); % get rid of empty waveform objects
    if numel(w)==0
        warning('no waveforms to plot')
        return
    end
    
    if ~exist('alignWaveforms', 'var')
            alignWaveforms = false;
    end
    
    if ~exist('subsetz','var')
        subsetz = [];
    end
    
    
    % get the first start time and last end time
    [starttimes endtimes]=gettimerange(w)
    snum = nanmin(starttimes);
    enum = nanmax(endtimes);
    
    
    % get the longest duration - in mode=='align' 
    durations = endtimes - starttimes;
    maxduration = nanmax(durations); 
    SECSPERDAY = 60 * 60 * 24;
    if numel(station) >1 || ~strcmp(station,'all')
        for ind=1:numel(station)
            staz = get(w,'station');
            st_ind(ind,:) = find(strcmp(station{ind},staz));
        end
        stzz=reshape(st_ind',1,numel(st_ind));
        w_subset=w(stzz);
        w=w_subset;
        max_vals=max_vals(stzz);
        min_vals=min_vals(stzz);
    end
    nwaveforms = numel(w);
    
	h=figure;
	trace_height=0.98/nwaveforms;
	left=0.125;
	width=0.8;
    %set(h, 'Position', [1000 1000 width+525 trace_height*nwaveforms+1000])
    set(h, 'Position', [1000 1000 width+1000 trace_height*nwaveforms+1400])
    for wavnum = 1:nwaveforms
        clear data
        clear dnum
        data=get(w(wavnum),'data');
        freq = get(w(wavnum), 'freq');
        dnum(1) = datenum(get(w(wavnum),'start'));
            for l = 2:numel(data)
                dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
            end 
        sta=get(w(wavnum),'station');
        chan=get(w(wavnum),'channel');
        ax(wavnum)=axes('Position',[left+0.025 0.95-wavnum*trace_height+0.05 width trace_height*.75]); 
        
        if alignWaveforms
            plot((dnum-min(dnum))*SECSPERDAY, data,'-k');
            set(gca, 'XLim', [0 maxduration*SECSPERDAY]);
        else
            plot(dnum, data,'-k'); 
            if strcmp(style, 'plotz')
                set(gca, 'XLim', [snum+datenum(0,0,0,0,0,4) enum-datenum(0,0,0,0,0,14)]); %10,8 PLSP %7,6 PL03
            else
                set(gca, 'XLim', [snum+datenum(0,0,0,0,0,3) enum-datenum(0,0,0,0,0,3)]);
            end
        end
        wavnum
        min_vals(wavnum)
        max_vals(wavnum)
        ylim([min_vals(wavnum) max_vals(wavnum)])
        ylabel(sprintf('%s\n%s ',sta,chan),'FontSize',9,'Rotation',90);

       
        
        %set(gca,'YTick',[],'YTickLabel',['']);
        datetick('x', 'keeplimits');
        if wavnum<nwaveforms
           set(gca,'XTickLabel',['']);
        end
        
       
        hold on
        indz = find(strcmp(oridStruct.(OrS).sta,sta));
        yl=ylim;
%         line([get(w(wavnum), 'EX_ARR_TIME'), get(w(wavnum), 'EX_ARR_TIME')], [yl(1), yl(2)], 'Color', 'k');
%         hold on
        if strcmp(style, 'linez')
            for val=1:numel(indz)
                line([oridStruct.(OrS).time_phase(indz(val)), oridStruct.(OrS).time_phase(indz(val))], [yl(1), yl(2)], 'Color', 'k', 'LineWidth', 2);
            end
            colors={'green','orange','purple','pink','silver'};
            %for david=1:5
            for david=1:numel(timez.(sta).refl.P)+1
                %if david <=4
                if david <=numel(timez.(sta).refl.P)
                    time_P = datenum(0,0,0,0,0,timez.(sta).refl.P(david))+min(dnum);
                    time_S = datenum(0,0,0,0,0,timez.(sta).refl.S(david))+min(dnum);
                    line([time_P, time_P], [yl(1), yl(2)], 'Color', rgb(colors{david}), 'LineWidth', 2);
                    line([time_S, time_S], [yl(1), yl(2)], 'Color', rgb(colors{david}), 'LineWidth', 2);
                else
                    line([datenum(0,0,0,0,0,timez.(sta).d.P)+min(dnum), datenum(0,0,0,0,0,timez.(sta).d.P)+min(dnum)], [yl(1), yl(2)], 'Color', rgb(colors{end}), 'LineWidth', 2);
                    line([datenum(0,0,0,0,0,timez.(sta).d.S)+min(dnum), datenum(0,0,0,0,0,timez.(sta).d.S)+min(dnum)], [yl(1), yl(2)], 'Color', rgb(colors{end}), 'LineWidth', 2);
                end
            end
        elseif strcmp(style, 'linez_shifted') % -----------------------in progress----------------------------------------
            for val=1:numel(indz)
                line([oridStruct.(OrS).time_phase(indz(val)), oridStruct.(OrS).time_phase(indz(val))], [yl(1), yl(2)], 'Color', 'k', 'LineWidth', 2);
            end
            try
                time_P = datenum(0,0,0,0,0,timez.(sta).d.P)+min(dnum);
                ind_P = intersect(find(strcmp(oridStruct.(OrS).sta,sta)),find(strcmp(oridStruct.(OrS).phase,'P')));
                time_Parr = oridStruct.(OrS).time_phase(ind_P);
                diffzP = time_Parr-time_P;
            end
            try
                time_S = datenum(0,0,0,0,0,timez.(sta).d.S)+min(dnum);
                ind_S = intersect(find(strcmp(oridStruct.(OrS).sta,sta)),find(strcmp(oridStruct.(OrS).phase,'S')))
                time_Sarr = oridStruct.(OrS).time_phase(ind_S);
                diffzS = time_Sarr-time_S;
            end
            colors={'green','orange','purple','pink','silver'};
            %for david=1:5
            for david=1:numel(timez.(sta).refl.P)+1
                %if david <=4
                if david <=numel(timez.(sta).refl.P)
                    %P
                    try
                    time_Pref = datenum(0,0,0,0,0,timez.(sta).refl.P(david))+min(dnum)
                    line([time_Pref+diffzP, time_Pref+diffzP], [yl(1), yl(2)], 'Color', rgb(colors{david}), 'LineWidth', 2);
                    %S
                    time_Sref = datenum(0,0,0,0,0,timez.(sta).refl.S(david))+min(dnum);
                    line([time_Sref+diffzS, time_Sref+diffzS], [yl(1), yl(2)], 'Color', rgb(colors{david}), 'LineWidth', 2);
                    end
                else
                    try
                    line([time_P+diffzP, time_P+diffzP], [yl(1), yl(2)], 'Color', rgb(colors{end}), 'LineWidth', 2);
                    line([time_S+diffzS, time_S+diffzS], [yl(1), yl(2)], 'Color', rgb(colors{end}), 'LineWidth', 2);
                    end
                end

                
            end
        elseif strcmp(style, 'plotz')
            for val=1:numel(indz)
                line([oridStruct.(OrS).time_phase(indz(val)), oridStruct.(OrS).time_phase(indz(val))], [yl(1), yl(2)], 'Color', 'k', 'LineWidth', 2);
            end
            data2=data(subsetz(1):subsetz(2));
            plot(dnum(subsetz(1):subsetz(2)), data2,'-r')
        end
        %data_construct = [dnum; data]'
%         if wavnum==1
%            title('','FontSize',10);
%         end
        %axis tight;
        
        % display mean on left, max on right
         a=axis;
%         tpos=a(1)+(a(2)-a(1))*.02; %horizontal position
         dpos=a(3)+(a(4)-a(3))*.7; %vertical position
%         text(tpos,dpos,sprintf('%5.0f',nanmean(data)),'FontSize',10,'Color','k');
% %         tpos=a(1)+(a(2)-a(1))*.4;
% %         text(tpos,dpos,sprintf(' %s',datestr(starttimes(wavnum),30)),'FontSize',10,'Color','k');
         tpos=a(1)+(a(2)-a(1))*.9;
         %text(tpos,dpos,sprintf('%5.0f',nanmax(abs(data))),'FontSize',10,'Color','k');

        
    end
%     if exist('ax','var')
%         linkaxes(ax,'x');
%         %samexaxis();
%         hlink = linkprop(ax,'XLim');
%         if ~alignWaveforms
%             %datetick('x', 'keeplimits');
%         end
%     end
xlabel('Time')

if numel(station)>1
    station={'subset'};
end
filename = sprintf('%s_%s %s_waveforms_%1.4f_%1.4f.png',OrS,style,station{1},fil(1),fil(2));
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');