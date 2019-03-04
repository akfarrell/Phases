% function phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(eq,stations_inEq{ctz},sec)
tic
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/reexternalre/')

addpath(genpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase'))
addpath('/Users/alexandrafarrell/Desktop/akfarrell/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Polarizemic-master/functions/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/fillPage/')
%clear;



close all; clc;
%fil=[2 25];
%fil=[10 20];
% [evidStruct, allevids] = get_eq_info();
load('siteStruct.mat');
%stylie = 'min'; %%%%% CHANGE!!!---------0----min,max,abz
warning('off')
%ch = {'HHE','HHN','HHZ'};
ch = {'HHZ','HHN','HHE'};
ch2 = {'HHZ','HHR','HHT'};
cutoff_val = 0.6;
S_pad = 25;
SECSPERDAY = 60 * 60 * 24;
num_samps = 31;
%close all
% eq = 2066;
% stations_inEq{ctz} = 'PLLL';
% eq = 2249;
% stations_inEq{ctz} = 'PLLL';
Pvel = 4.2; %km/s
Svel = 2.35; %km/s
%
eqs = [1811,1816, 1825, 1830, 1831, 1832,1833, 1834, 1835, 1836, 1839, 1847, 1848,...
    1855, 1859, 1861, 1862, 1864, 1866, 1868, 1870, 1871, 1872, 1874, 1876,...
    1878, 1881, 1882, 1883, 1943, 1944, 1945, 1954, 1955, 1957, 1968, 1969, ...
    1970, 1974, 1987, 1988, 1987, 1988, 1992, 1993, 1994, 1996, 1997, 1998, ...
    1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,...
    2011, 2012, 2013, 2014, 2015, 2018, 2019, 2020, 2021, 2022, 2023, 2023, ...
    2024, 2025, 2026, 2027, 2031, 2033, 2035, 2050, 2051, 2055, 2058, 2060, ...
    2061, 2068, 2080, 2089, 2090, 2097, 2110, 2112, 2116, 2120, 2132, 2133, ...
    2136, 2142, 2156, 2157, 2160, 2165, 2166, 2167, 2171, 2177, 2207, 2215, ...
    2233, 2244, 2249, 2256, 2264, 2267, 2284, 2285, 2293, 2297, 2298, 2300, ...
    2301, 2305, 2306, 2307, 2312, 2313, 2314, 2321, 2322, 2338, 2340, 2360, ...
    2364, 2368, 2369, 2370, 2372, 2376, 2393, 2403, 2404, 2405, 2409, 2411, ...
    2414, 2429, 2432, 2433, 2436, 2437, 2441, 2442, 2458, 2463, 2468, 2470,...
    2482, 2484, 2488, 2489, 2490, 2492, 2493, 2499, 2501, 2543, 2562, 2585];

for ctz1 = 1:numel(eqs)
    % eq = 2501;
    % stations_inEq{ctz} = 'PLMN';
    eq = eqs(ctz1);
    [evidStruct, allevids] = get_eq_info();
    
    Pvel = 4.2; %km/s
    Svel = 2.35; %km/s
    %val_try = 569; % If commented out, find val_try (is P_ind or S_ind) ~line
    %110
    count=find(allevids == eq);
    if isempty(count)
        clear evidStruct allevids count
        [evidStruct, allevids] = get_eq_info('/Users/alexandrafarrell/Desktop/dbplutons_alex','true'); %not relocated, with errors
        count = find(allevids == eq);
    end
    if isempty(count)
        continue
    end
    
    %create and clean waveform object
    
    
    [w_raw,EvS,stations_inEq] = get_wf(allevids(count),evidStruct);
    %     if isempty(find(strcmp(get(w_raw,'channel'),'HHE')))
    %         error('Stations not aligned correctly')
    %     end
    %w_clean_filt = waveform_clean(w_raw, filterobject('b', fil, 2));
    w_clean = waveform_clean(w_raw);
    %save(filename_wPath,'w_clean', 'EvS', 'stations_inEq');
    %end
    stationz = get(w_clean,'station');
    
    for ctz = 1:numel(stations_inEq)
        count2 = min(find(strcmp(stations_inEq{ctz},stationz)));%:3:numel(w_clean) %43-45 is PLMN 43:3:43
        chan_temp = get(w_clean,'channel');
        if isempty(find(strcmp(chan_temp(count2:count2+2),'HHE')))
            fprintf('Stations not aligned correctly %d', eq)
            continue
        end
        %HHE = count, HHN = count+1, HHZ = count+2
        clear P_ind ind_P time_Parr ind_S time_Sarr
        ind_P = intersect(find(strcmp(evidStruct.(EvS).sta,stationz{count2})),find(strcmp(evidStruct.(EvS).phase,'P')));
        if isempty(ind_P)
            fprintf('Station does not have P wave %s %d', stations_inEq{ctz},eq)
            continue
        end
        
        time_Parr = evidStruct.(EvS).time_phase(ind_P);
        
        try
            ind_S = intersect(find(strcmp(evidStruct.(EvS).sta,stationz{count2})),find(strcmp(evidStruct.(EvS).phase,'S')));
            time_Sarr = evidStruct.(EvS).time_phase(ind_S);
        catch
            ind_S = numel(get(w_clean(1),'data'));
        end
        dnum = zeros(1,numel(get(w_clean(count2),'data')));
        dnum(1) = datenum(get(w_clean(count2),'start'));
        freq = get(w_clean(count2), 'freq');
        for l = 2:numel(get(w_clean(count2),'data'))
            dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
        end
        P_ind = find(dnum>=time_Parr, 1,'first');
        if isempty(P_ind)
            P_ind = 0;
        end
        S_ind = find(dnum>=time_Sarr, 1,'first');
        if isempty(S_ind)
            S_ind = numel(get(w_clean(1),'data'))-num_samps;
        end
        %%
        % fname = sprintf('corr_%d_%s_%s.mat',allevids(count),stationz{count2},stylie{1});
        % directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/wf_objs';
        % directory2 = sprintf('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/corrs/%d',allevids(count));
        % filename_wPath2 = fullfile(directory2,fname);
        % load(filename_wPath2)
        % c_orig=c;
        %P_pad = 15
        
        
        numz = find(strcmp(siteStruct.sta,stations_inEq{ctz}));
        % find azimuth, distance, and incidence angle
        az = azimuth(siteStruct.lat(numz), siteStruct.lon(numz),evidStruct.(sprintf('eq_%d',eq)).lat(1),evidStruct.(sprintf('eq_%d',eq)).lon(1));
        starttime = evidStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(evidStruct.(sprintf('eq_%d',eq)).phase,'P')),...
            find(strcmp(evidStruct.(sprintf('eq_%d',eq)).sta,stations_inEq{ctz}))))-datenum(0,0,0,0,0,1);
        dist = distance(siteStruct.lat(numz), siteStruct.lon(numz),evidStruct.(sprintf('eq_%d',eq)).lat(1),evidStruct.(sprintf('eq_%d',eq)).lon(1))*111.12;
        inc = atand((siteStruct.elev(numz)+evidStruct.(sprintf('eq_%d',eq)).depth(1))/dist);
        
        az2 = az;
        if az > 180
            az = az-180;
        else
            az = az+180;
        end
        
        
        
        val_try = P_ind; indie = 1; ph_used = 'P';
        
        
        endtime = P_ind+ceil((dist/Svel-dist/Pvel)*100)+100;
        % Find best filter using periodogram and filter to that
        fs = get(w_clean(count2),'freq');
        w_clean_copy = w_clean;
        w_clean_copy = w_clean_copy(count2:count2+2);
        compnt = get(w_clean_copy,'channel');
        data = get(w_clean_copy,'data');
        %close all
        errorz = zeros(3,2);
        for countz = 1:numel(data)
            sig_data.(compnt{countz}) = data{countz}(val_try-5:val_try+36);
            %added 5 samples before val_try to allow for a padding for fft
            %also added 5 samples after 32 samples, to allow for padding for fft
            leng = numel(sig_data.(compnt{countz}));
            nfft = leng;
            %[pdgrm2.(compnt{countz}),f_pdgrm2.(compnt{countz})] = periodogram(sig_data.(compnt{countz}),[],nfft,fs);
            [pdgrm1.(compnt{countz}),f_pdgrm1.(compnt{countz})] = periodogram(sig_data.(compnt{countz}),hamming(leng),nfft,fs);
            [pdgrm.(compnt{countz}),f_pdgrm.(compnt{countz})] = periodogram(detrend(sig_data.(compnt{countz})),hamming(leng),nfft,fs);
            
            %CHECK WITH NORMALIZING AND ADDING TOGETHER SPECTRA
            max_p = max(pdgrm.(compnt{countz}));
            pdgrm_norm.(compnt{countz}) = pdgrm.(compnt{countz})./max_p;
            clear max_p
            
            % TRY THIS WITH DETERMINE PEAKS OF BOTH!!!!!
            try
                [k.(compnt{countz}),w.(compnt{countz})] = det_pks(pdgrm.(compnt{countz}),f_pdgrm.(compnt{countz}),3,'pdgrm','low-freq')
            catch
                [k.(compnt{countz}),w.(compnt{countz})] = det_pks(pdgrm.(compnt{countz}),f_pdgrm.(compnt{countz}),3,'pdgrm')
                disp('NOT LOW-FREQ')
            end
            %[k1.(compnt{countz}),w1.(compnt{countz})] = det_pks(pdgrm1.(compnt{countz}),f_pdgrm1.(compnt{countz}),3,'pdgrm')
            %w1.(compnt{countz})=w1.(compnt{countz})/2;
            w.(compnt{countz})=w.(compnt{countz})/2;
            errorz(countz,:) = [k.(compnt{countz})-w.(compnt{countz}),k.(compnt{countz})+w.(compnt{countz})];
            %errorz1.(compnt{countz}) = [k1.(compnt{countz})+w1.(compnt{countz}),k1.(compnt{countz})-w1.(compnt{countz})];
            
        end
        % Determine peak frequency of normalized, composite spectra
        
        pdgrm_tot = pdgrm_norm.HHZ+pdgrm_norm.HHE+pdgrm_norm.HHN;
        [pk_tot,w_tot,which_pk] = det_pks(pdgrm_tot,f_pdgrm.HHZ,3,'pdgrm','low-freq')
        if strcmp(stations_inEq{ctz},'PLCM') && eq==2165 && indie == 1
            warning('on','all')
            warning('PLCM for eq 2165 chooses 40 Hz - forcing to 14.29')
            pk_tot = 14.29; %PLCM picks dominant frequency at 40 Hz, because the first arrival is shitty
            which_pk = '-';
            warning('off','all')
        end
        
        %% Check into values for best fitting frequencies and see how they compare to waveforms
        intersect_vals = RangesIntersect(errorz(1,:), errorz(2,:), errorz(3,:));
        a = linspace(errorz(1,1),errorz(1,2),50);
        b = linspace(errorz(2,1),errorz(2,2),50);
        c = linspace(errorz(3,1),errorz(3,2),50);
        d = mean(horzcat(a,b,c));
        
        %make sine waves of frequency picked up by det_pks
        sig_len = numel(sig_data.HHZ);
        x = 0:sig_len-1; %128(127) okay, 256(255) is best
        fs = 100;
        y = 100*sin(2*pi*d/fs*x);
        y2 = 100*sin(2*pi*pk_tot/fs*x);
        %         figure() %figure to show what the proposed frequency looks like comared to the waveforms
        %         plot(x,y)
        %         hold on
        %         plot(x,y2*1000,'b--')
        %         plot(x,sig_data.HHZ,'r')
        %         plot(x,sig_data.HHN,'k')
        %         plot(x,sig_data.HHE,'g')
        %         title(sprintf('d - pkTot = %2.2f',d-pk_tot))
        %         hold off
        %         %%
        %         fprintf('d - pkTot = %2.2f\n',d-pk_tot)
        %         disp('using pk_tot')
        
        %% Use the filter above on the data
        fil = [0.8 min(max(25,pk_tot+10),40)];
        % fil = [max(pk_tot-10,0.5) pk_tot+10];%[max(pk_tot-5,0.5) pk_tot+5];%[0.5 25];%[pk_tot-2 pk_tot+2];
        w_clean_filt = waveform_clean(w_clean, filterobject('b', fil, 2));
        w_clean_nofilt = w_clean;
        % warning('on','all')
        % warning('USING UNFILTERED DATA IN W_CLEAN_FILT')
        % warning('off','all')
        % w_clean_filt = w_clean_nofilt;
        % %This shows the ellipticity, azimuth, and incidence of the ray at that
        %filter
        
        %% Subset data to around the P- and S-wave picks
        
        % try
        %     endtime = evidStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(evidStruct.(sprintf('eq_%d',eq)).phase,'S')),...
        %         find(strcmp(evidStruct.(sprintf('eq_%d',eq)).sta,stations_inEq{ctz}))))+datenum(0,0,0,0,0,1);
        %     w_clean_e = extract(w_clean_filt,'TIME',starttime,endtime);
        % catch
        %     w_clean_e = extract(w_clean_filt,'TIME',starttime);
        % end
        
        
        
        
        
        %%
        for count3 = 1:numel(ch)
            %W(1,3-count3+1) = w_clean_e(count2+count3-1);
            W_all(1,3-count3+1) = w_clean_filt(count2+count3-1);
            %organizing and finding only waveforms for event of interest
            W_all_unfilt(1,3-count3+1) = w_clean_nofilt(count2+count3-1);
            
        end
        %         try
        %             W_all_short = extract(W_all,'TIME', time_Parr-datenum(0,0,0,0,0,2), time_Sarr+datenum(0,0,0,0,0,2));
        %             W_all_unfilt_short = extract(W_all_unfilt,'TIME', time_Parr-datenum(0,0,0,0,0,2), time_Sarr+datenum(0,0,0,0,0,2));
        %         catch
        %             W_all_short = extract(W_all,'TIME', time_Parr-datenum(0,0,0,0,0,3), dnum(endtime));
        %             W_all_unfilt_short = extract(W_all_unfilt,'TIME', time_Parr-datenum(0,0,0,0,0,3), dnum(endtime));
        %         end
        clear count3
        %% rotate waveforms
        wndo = 32;%30;%20;%14;%31;
        TC=threecomp(W_all,az2);
        %         TCS = threecomp(W_all_short,az2);
        TCR = TC.rotate();
        %         TCSR = TCS.rotate();
        TCW = TCR.waveform();
%         TCRW = TCSR.waveform();
        
        
%         TCUS = threecomp(W_all_unfilt_short,az2); %az2
%         TCUSR = TCUS.rotate();
        
        
        
        
        
        wind_use = max(wndo,ceil(200/pk_tot));
        fprintf('Using window length %2.4f\n', wind_use);
        [TCP,vals,vecs] = TCR.particlemotion(0.01,wind_use/100);
        [TCP_overall,vals_overall,vecs_overall] = TCR.particlemotion();
        %[TCP,vals,vecs] = TCR.particlemotion();
        [lambda, I]=sort(vals,1,'descend');
        X = vecs(:,I(1),:); %get greatest eigenvector
        %         TCSRP = TCSR.particlemotion(0.01,wind_use/100);
        %         TCSRPL = TCSR.particlemotion();
%         TCUSRP = TCUSR.particlemotion(0.01,wind_use/100);
%         TCUSRPL = TCUSR.particlemotion();
        
        TC_unfilt = threecomp(W_all_unfilt,az);
        %TCR_unfilt = TC_unfilt.rotate();
        TCW_unfilt = TC_unfilt.waveform();
        [TCP_unfilt, vals_unfilt, vecs_unfilt] = TC_unfilt.particlemotion(0.01,wind_use/100);
        tcrp_az = get(TCP,'azimuth');
        tcrp_inc = get(TCP,'inclination');
        
        tcrp_az_u = get(TCP_unfilt,'azimuth');
        tcrp_inc_u = get(TCP_unfilt,'inclination');
        %phase_plots(Haystack_data, P_ind, S_ind, get(tcrp_az,'data'),get(tcrp_inc,'data'),get(tcrp_rec,'data'),eq,stasz);
        
        tcrp_az_data = get(tcrp_az,'data');
        if numel(tcrp_az_data) ~= numel(data{1})
            continue
        end
        
        %% Check and plot directionality and polarity of waveforms
        %         dtac_i = [Haystack_data.(ch{1})'; Haystack_data.(ch{3})'; Haystack_data.(ch{2})'];
        %         dtac_unfilt = [Haystack_data_unfilt.(ch{1})'; Haystack_data_unfilt.(ch{2})'; Haystack_data_unfilt.(ch{3})'];
        %         dtac = dtac_i(:,max(P_ind-100,1):S_ind+100);
        % dtac = dtac_i;
        sps = 100;
        % window length in (s)
        windt = wndo/100;
        % fraction of window for subwindows used in cross-spectrum estimation
        fsw = 2/wndo;
        % fraction of window for overlap between time bins
        fov = 0.9;
        
        % time-frequency polarization analysis
        %         [azim, incd, ellip, rpol, ppol, tvecC_zs, F, Cxypzz, Cxypnn, Cxypee, largest_eig] = ...
        %             polartf_cross_spectrum(dtac,windt,sps,fsw,fov);
        %         largest_eig;
        
        %
        %         [aa bb] = size(Cxypzz);
        %         pctr = 0; %0.002;
        %         %mxc = max(max(Cxypzz));
        %         mxc = max(max([Cxypzz Cxypee Cxypnn]));
        %         for ii=1:aa
        %             for jj=1:bb
        %                 if ((Cxypzz(ii,jj)> pctr*mxc) | (Cxypnn(ii,jj)> pctr*mxc) | (Cxypee(ii,jj)> pctr*mxc))
        %                     incdm(ii,jj) = incd(ii,jj);
        %                     azimm(ii,jj) = azim(ii,jj);
        %                     ellipm(ii,jj) = ellip(ii,jj);
        %                 else
        %                     incdm(ii,jj) = NaN;
        %                     azimm(ii,jj) = NaN;
        %                     ellipm(ii,jj) = NaN;
        %                 end
        %             end
        %         end
        
        %         rec_vals = get(tcrp_rec,'data');
        %         plan_vals = get(tcrp_plan,'data');
        %         [valrec,indrec] = max(rec_vals(val_try:val_try+num_samps));
        %         if plan_vals(indrec+val_try-1)<0.5
        %             disp('Low degree of planarity!!')
        %         end
        ph_indz = 1;
%         Haystack_data.HHE = get(w_clean_filt(count2),'data');
%         Haystack_data.HHN = get(w_clean_filt(count2+1),'data');
%         Haystack_data.HHZ = get(w_clean_filt(count2+2),'data');
%         tcrp_rec = zeros(size(tcrp_az));
%         tcrp_plan = tcrp_rec;
%         tcrp_en = tcrp_rec;
%         plot_threecomp_directionality(endtime,val_try,Haystack_data,P_ind,S_ind,...
%             tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az,inc,tcrp_az_data,ph_indz)
        tcrp_inc_data = get(tcrp_inc,'data');
        tcrp_inc_data = 90-tcrp_inc_data;
        
        clear orientation_struct
        load('orientation_struct.mat')
        orientation_struct.(sprintf('eq_%d',eq)).(stations_inEq{ctz}).dist = dist;
        orientation_struct.(sprintf('eq_%d',eq)).(stations_inEq{ctz}).az = az;
        orientation_struct.(sprintf('eq_%d',eq)).(stations_inEq{ctz}).az2 = az2;
        orientation_struct.(sprintf('eq_%d',eq)).(stations_inEq{ctz}).incWest = mean(tcrp_inc_data(val_try:val_try+6));
        orientation_struct.(sprintf('eq_%d',eq)).(stations_inEq{ctz}).azWest = mean(tcrp_az_data(val_try:val_try+6));
        orientation_struct.(sprintf('eq_%d',eq)).(stations_inEq{ctz}).std_incWest = std(tcrp_inc_data(val_try:val_try+6));
        orientation_struct.(sprintf('eq_%d',eq)).(stations_inEq{ctz}).std_azWest = std(tcrp_az_data(val_try:val_try+6));
        orientation_struct.(sprintf('eq_%d',eq)).(stations_inEq{ctz}).freq = pk_tot;
        save('orientation_struct.mat','orientation_struct')
    end
    clearvars -except eqs ctz1 Pvel Svel num_samps SECSPERDAY S_pad cutoff_val ch...
        ch2 siteStruct
    close all
end
toc