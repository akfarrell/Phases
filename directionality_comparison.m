% determine the differences between the expected azimuths and the
% calculated ones, add that to phaseStruct
load('phaseStruct.mat')
events = fieldnames(phaseStruct);
for count = 1:numel(events)
    stations = fieldnames(phaseStruct.(events{count}));
    for count2 = 1:numel(stations)
        west_az = phaseStruct.(events{count}).(stations{count2}).azWest(1);
        haney_az = phaseStruct.(events{count}).(stations{count2}).haney_az(1);
        az = phaseStruct.(events{count}).(stations{count2}).azRange;
        if haney_az <= 90
            haney_az = 90-haney_az;
        elseif haney_az <= 180
            haney_az = 270-haney_az;
        end
        
        if west_az <=180
            west_az(2) = west_az+180;
            west_az(3) = west_az(1)+360;
        else
            west_az(2) = west_az-180;
            west_az(3) = west_az(1)+180;
            
        end
        
        if haney_az <=180
            haney_az(2) = haney_az+180;
            haney_az(3) = haney_az(1)+360;
        else
            haney_az(2) = haney_az-180;
            haney_az(3) = haney_az(1)+180;
        end
        west_az;
        haney_az;
%         try
            haney_diffs = [haney_az(3)-az(1), haney_az(3)-az(2), haney_az(1)-az(1), haney_az(1)-az(2), haney_az(2)-az(1), haney_az(2)-az(2)];
            west_diffs = [west_az(3)-az(1), west_az(3)-az(2), west_az(1)-az(1), west_az(1)-az(2), west_az(2)-az(1), west_az(2)-az(2)];
%         catch
%             haney_diffs = [haney_az(1)-az(1), haney_az(1)-az(2), haney_az(2)-az(1), haney_az(2)-az(2)];
%             west_diffs = [west_az(1)-az(1), west_az(1)-az(2), west_az(2)-az(1), west_az(2)-az(2)];
%         end
        [val,ind] = min(abs(haney_diffs));
        haney_final_diff = haney_diffs(ind);
        clear val ind
        [val,ind] = min(abs(west_diffs));
        west_final_diff = west_diffs(ind);
        if abs(west_final_diff) >170
            west_final_diff
            west_diffs
            west_az
            az
            error('west solution sucks')
        end
        if abs(haney_final_diff) > 170
            haney_final_diff
            haney_diffs
            haney_az
            az
            error('haney solution sucks')
        end
        
        phaseStruct.(events{count}).(stations{count2}).westVals = west_az;
        phaseStruct.(events{count}).(stations{count2}).haneyVals = haney_az;
        phaseStruct.(events{count}).(stations{count2}).haney_final_diff = haney_final_diff;
        phaseStruct.(events{count}).(stations{count2}).west_final_diff = west_final_diff;
    end
    clearvars -except count phaseStruct events
end
save('phaseStruct.mat','phaseStruct')