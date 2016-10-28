function phase_mulplt(w, alignWaveforms, max_vals, min_vals)
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
    
    if ~exist('id', 'var')
        id = 'eq';
    end
    
    % get the first start time and last end time
    [starttimes endtimes]=gettimerange(w);
    snum = nanmin(starttimes);
    enum = nanmax(endtimes);
    
    
    % get the longest duration - in mode=='align' 
    durations = endtimes - starttimes;
    maxduration = nanmax(durations); 
    SECSPERDAY = 60 * 60 * 24;
    
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
            if strcmp(id, 'selected')
                if strcmp(chan,'HHZ') == 1 || strcmp(chan, 'HHE') == 1
                    plot(dnum, data,'-k');
                end
            else
              plot(dnum, data,'-k');  
            end
            set(gca, 'XLim', [snum enum]);
        end
        ylim([min_vals(wavnum) max_vals(wavnum)])
        if strcmp(id,'selected')
            if strcmp(chan,'HHZ') == 1 || strcmp(chan, 'HHE') == 1
                if strcmp(sta,'PL07')==1 || strcmp(sta,'PLSE') == 1 || strcmp(sta,'PLMK') ==1
                   ylabel(sprintf('%s\n%s ',sta,chan),'FontSize',9,'Rotation',90, 'FontWeight', 'bold');
                else
                   ylabel(sprintf('%s\n%s ',sta,chan),'FontSize',9,'Rotation',90); 
                end
            end
        else
            ylabel(sprintf('%s\n%s ',sta,chan),'FontSize',9,'Rotation',90);
        end
       
        
        %set(gca,'YTick',[],'YTickLabel',['']);
        datetick('x', 'keeplimits');
        if wavnum<nwaveforms
           set(gca,'XTickLabel',['']);
        end
        
       
%         hold on
%         yl=ylim;
%         line([get(w(wavnum), 'EX_ARR_TIME'), get(w(wavnum), 'EX_ARR_TIME')], [yl(1), yl(2)], 'Color', 'k');
%         hold on
        %line([time_value, time_value], [yl(1), yl(2)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
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
xlabel('Time');

%filename = sprintf('%s_%s_waveforms_%1.4f_%1.4f.png',name,id,fil(1),fil(2));
%hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');