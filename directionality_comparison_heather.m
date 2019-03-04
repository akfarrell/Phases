% determine the differences between the expected azimuths and the
% calculated ones, add that to orientation_struct
load('orientation_struct.mat')
events = fieldnames(orientation_struct);
for count = 1:numel(events)
    stations = fieldnames(orientation_struct.(events{count}));
    for count2 = 1:numel(stations)
        west_az = orientation_struct.(events{count}).(stations{count2}).azWest(1);
        az = orientation_struct.(events{count}).(stations{count2}).azRange;
        
        if west_az <=180
            west_az(2) = west_az+180;
            west_az(3) = west_az(1)+360;
            west_az(4) = west_az(1)-180;
        else
            west_az(2) = west_az-180;
            west_az(3) = west_az(1)+180;
            west_az(4) = west_az(2)-180;
            
        end
        
        west_az;
        %         try
        west_diffs = [west_az(3)-az(1), west_az(3)-az(2), west_az(1)-az(1), west_az(1)-az(2), west_az(2)-az(1), west_az(2)-az(2),...
            west_az(4)-az(1), west_az(4)-az(2)];
        %         catch
        %             haney_diffs = [haney_az(1)-az(1), haney_az(1)-az(2), haney_az(2)-az(1), haney_az(2)-az(2)];
        %             west_diffs = [west_az(1)-az(1), west_az(1)-az(2), west_az(2)-az(1), west_az(2)-az(2)];
        %         end
        
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
        
        orientation_struct.(events{count}).(stations{count2}).westVals = west_az;
        orientation_struct.(events{count}).(stations{count2}).west_final_diff = west_final_diff;
    end
    clearvars -except count orientation_struct events
end
save('orientation_struct.mat','orientation_struct')