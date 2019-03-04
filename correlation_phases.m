clear; close all; clc;
addpath('/Users/alexandrafarrell/Desktop/akfarrell/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Polarizemic-master/functions/')
[evidStruct, allevids] = get_eq_info('/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged'); %just relocated, no errors
[evidStruct_error, allevids_error] = get_eq_info('/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged','true'); %relocated, with errors (excludes all m<0.5

[evidStruct_unchecked_error, allevids_unchecked_error] = get_eq_info('/Users/alexandrafarrell/Desktop/dbplutons_alex','true'); %not relocated, with errors
load('phaseStruct.mat')
%%
% evid_test = [2120, 2121, 2122, 2123, 2124, 2125, 2126, 2127, 2128, 2130, 2131, 2132];
% evid_test = [1992 1993 1994 1995 1996 1997];
evid_test = [1847 1848 1849 1850 1851 1852 1853 1854 1855 1856 1857 1858 1859 1860 1861 1862 1863 1864 1866 1867 1868 1869 1870 1871 1872 1874 1875];
% evid_test = [2008 2009 2010 2011 2012 2013 2014 2015 2016 2018 2019 2020 2021 2022 2023 2024 2025 2026];
% evid_test = [2305 2306 2307 2308 2309 2310 2311 2312];
% evid_test = [2409 2410 2411 2412 2413];

time_num = zeros(1,numel(evid_test));
time_num_PLLL = time_num;
sta = {'PLSS'};
sta2 = {'PLLL'};
for count = 1:numel(evid_test)
    ind_P = intersect(find(strcmp(evidStruct.(sprintf('eq_%d',evid_test(count))).sta,sta)),find(strcmp(evidStruct.(sprintf('eq_%d',evid_test(count))).phase,'P')));
%     ind_P = find(strcmp(evidStruct.(sprintf('eq_%d',evid_test(count))).phase,'P'));
    time_num(count) = evidStruct.(sprintf('eq_%d',evid_test(count))).time_phase(ind_P);
    %clear ind_P
    ind_P_PLLL = intersect(find(strcmp(evidStruct.(sprintf('eq_%d',evid_test(count))).sta,sta2)),find(strcmp(evidStruct.(sprintf('eq_%d',evid_test(count))).phase,'P')));
%     ind_P = find(strcmp(evidStruct.(sprintf('eq_%d',evid_test(count))).phase,'P'));
    time_num_PLLL(count) = evidStruct.(sprintf('eq_%d',evid_test(count))).time_phase(ind_P_PLLL);
    
end

ds = datasource('antelope','/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged');

scnl_PLSS = scnlobject(sta,'HHZ');
c_PLSS = correlation(ds,scnl_PLSS,time_num,-3,5)

scnl_PLLL = scnlobject(sta2,'HHZ');
c_PLLL = correlation(ds,scnl_PLLL,time_num_PLLL,-3,5)

%%
c_PLLL = taper(c_PLLL);
c_PLSS = taper(c_PLSS);
c_PLLL = butter(c_PLLL, [0.8 26.6]);

c_PLLL =xcorr(c_PLLL,[0 2]);
c_PLSS = xcorr(c_PLSS, [0 1]);

%%
c_PLLL = sort(c_PLLL);
plot(c_PLLL,'corr');
set(gcf,'Position',[50 50 500 400]);
corr_matrix_PLLL = get(c_PLLL,'CORR');

c_PLLL = massage_xcorr(c_PLLL);
c_PLSS = massage_xcorr(c_PLSS);

%%
PLLL_clust = get(c_PLLL,'CLUST');
PLSS_clust = get(c_PLSS,'CLUST');
vals = unique(PLLL_clust);
for count2 = 1:numel(unique(PLLL_clust))
    ind = find(PLLL_clust==vals(count2))
    evid_subset_PLLL.(sprintf('cluster_%d',count2)) = evid_test(ind)
    clear ind
end
clear vals count2

vals = unique(PLSS_clust);
for count2 = 1:numel(unique(PLSS_clust))
    ind = find(PLSS_clust==vals(count2))
    evid_subset_PLSS.(sprintf('cluster_%d',count2)) = evid_test(ind)
    clear ind
end
clear vals count2

