function makeTXTfiles(phaseStruct,eq,stations,testvar)
    if strcmp(testvar,'az_W_calc')
        fid1 = fopen(sprintf('%s_diff_expP_auto.txt',eq),'w');
        fprintf(fid1,'Vals-expected P-wave arrival autocalc\ndiff_expP_auto.txt from maps_phasestuff.m\n\n');
        
        fid2 = fopen(sprintf('%s_diff_calcP_auto.txt',eq),'w');
        fprintf(fid2,'Vals-calculated P-wave arrival autocalc\ndiff_calcP_auto.txt from maps_phasestuff.m\n\n');
        
        fid3 = fopen(sprintf('%s_just_vals_auto.txt',eq),'w');
        fprintf(fid3,'Raw values autocalc\njust_vals_auto.txt from maps_phasestuff.m\n\n');
    elseif strcmp(testvar,'azWest')
        fid1 = fopen(sprintf('%s_diff_expP.txt',eq),'w');
        fprintf(fid1,'Vals-expected P-wave arrival\ndiff_expP.txt from maps_phasestuff.m\n\n');
        
        fid2 = fopen(sprintf('%s_diff_calcP.txt',eq),'w');
        fprintf(fid2,'Vals-calculated P-wave arrival\ndiff_calcP.txt from maps_phasestuff.m\n\n');
        
        fid3 = fopen(sprintf('%s_just_vals.txt',eq),'w');
        fprintf(fid3,'Raw values\njust_vals.txt from maps_phasestuff.m\n\n');
    elseif strcmp(testvar,'inc_W_calc')
        fid1 = fopen(sprintf('%s_diff_inc_expP_auto.txt',eq),'w');
        fprintf(fid1,'Vals-expected P-wave arrival inc autocalc\ndiff_inc_expP_auto.txt from maps_phasestuff.m\n\n');
        
        fid2 = fopen(sprintf('%s_diff_inc_calcP_auto.txt',eq),'w');
        fprintf(fid2,'Vals-calculated P-wave arrival inc autocalc\ndiff_inc_calcP_auto.txt from maps_phasestuff.m\n\n');
        
        fid3 = fopen(sprintf('%s_just_inc_vals_auto.txt',eq),'w');
        fprintf(fid3,'Raw values inc autocalc\njust_inc_vals_auto.txt from maps_phasestuff.m\n\n');
    elseif strcmp(testvar,'incWest')
        fid1 = fopen(sprintf('%s_diff_inc_expP.txt',eq),'w');
        fprintf(fid1,'Vals-expected P-wave arrival inc\ndiff_inc_expP.txt from maps_phasestuff.m\n\n');
        
        fid2 = fopen(sprintf('%s_diff_inc_calcP.txt',eq),'w');
        fprintf(fid2,'Vals-calculated P-wave arrival inc\ndiff_inc_calcP.txt from maps_phasestuff.m\n\n');
        
        fid3 = fopen(sprintf('%s_just_inc_vals.txt',eq),'w');
        fprintf(fid3,'Raw values inc\njust_inc_vals.txt from maps_phasestuff.m\n\n');
    end
    
    fprintf(fid1,'Stn\tExp\tP\tP1\tP2\tP3\tS');
    fprintf(fid2,'Stn\tExp\tP\tP1\tP2\tP3\tS');
    fprintf(fid3,'Stn\tExp\tP\tP1\tP2\tP3\tS');
    for count = 1:numel(stations)
        vals = phaseStruct.(eq).(stations{count}).(testvar);
        if strfind(testvar,'az')
            if phaseStruct.(eq).(stations{count}).az_exp_adj > 360
                expected = phaseStruct.(eq).(stations{count}).az_exp_adj-360;
            else
                expected = phaseStruct.(eq).(stations{count}).az_exp_adj;
            end
        elseif strfind(testvar,'inc')
            expected = phaseStruct.(eq).(stations{count}).inc_exp;
        end
        
        for count2 = 1:numel(vals)
            if vals(count2) > 360
                vals(count2) = vals(count2)-360;
            end
        end
        if isequal(phaseStruct.(eq).(stations{count}).phase,{'P','P1','S'})
            fprintf(fid1,'\n%s\t%3.1f\t%3.1f\t%3.1f\t-\t-\t%3.1f',stations{count},expected,vals(1)-expected,vals(2)-expected,vals(3)-expected); %exp
            fprintf(fid2,'\n%s\t%3.1f\t%3.1f\t%3.1f\t-\t-\t%3.1f',stations{count},expected-vals(1),vals(1),vals(2)-vals(1),vals(3)-vals(1)); %calc
            fprintf(fid3,'\n%s\t%3.1f\t%3.1f\t%3.1f\t-\t-\t%3.1f',stations{count},expected,vals(1),vals(2),vals(3)); %just vals
        elseif isequal(phaseStruct.(eq).(stations{count}).phase,{'P','P1'})
            fprintf(fid1,'\n%s\t%3.1f\t%3.1f\t%3.1f\t-\t-\t-',stations{count},expected,vals(1)-expected,vals(2)-expected); %exp
            fprintf(fid2,'\n%s\t%3.1f\t%3.1f\t%3.1f\t-\t-\t-',stations{count},expected-vals(1),vals(1),vals(2)-vals(1)); %calc
            fprintf(fid3,'\n%s\t%3.1f\t%3.1f\t%3.1f\t-\t-\t-',stations{count},expected,vals(1),vals(2)); %just vals
        elseif isequal(phaseStruct.(eq).(stations{count}).phase,{'P','P1','P2','P3','S'})
            fprintf(fid1,'\n%s\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f',stations{count},expected,vals(1)-expected,vals(2)-expected,vals(3)-expected, vals(4)-expected, vals(5)-expected); %exp
            fprintf(fid2,'\n%s\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f',stations{count},expected-vals(1),vals(1),vals(2)-vals(1),vals(3)-vals(1),vals(4)-vals(1),vals(5)-vals(1)); %calc
            fprintf(fid3,'\n%s\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f',stations{count},expected,vals(1),vals(2),vals(3),vals(4),vals(5)); %just vals
        elseif isequal(phaseStruct.(eq).(stations{count}).phase,{'P','P1','P2','S'})
            fprintf(fid1,'\n%s\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t-\t%3.1f',stations{count},expected,vals(1)-expected,vals(2)-expected,vals(3)-expected,vals(4)-expected); %exp
            fprintf(fid2,'\n%s\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t-\t%3.1f',stations{count},expected-vals(1),vals(1),vals(2)-vals(1),vals(3)-vals(1),vals(4)-vals(1)); %calc
            fprintf(fid3,'\n%s\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t-\t%3.1f',stations{count},expected,vals(1),vals(2),vals(3),vals(4)); %just vals
        else
            warning('This phase combo is not represented')
        end
        clear vals
    end
    fclose(fid1);
    fclose(fid2);
    fclose(fid3);
end