function siteStruct = siteStruct_for_Heather()
databasePath='/raid/data/antelope/databases/PLUTONS/dbmerged';



%siteStruct = struct('sta', [], 'ondate', [], 'offdate', [], 'lat', [], 'lon', [], 'elev', []);
%Create siteStruct structure to be filled in later   


dbp = dbopen(databasePath, 'r'); %Open to read database at specified databasePath
db = dblookup_table(dbp, 'site'); %Look up origin table, assign that to variable db


[sta, ondate, offdate, lat, lon, elev] = ...
dbgetv(db, 'sta','ondate','offdate','lat','lon','elev'); %assign variables to each column specified    
dbclose(db);
siteStruct = struct('sta', {sta},'ondate',ondate,'offdate',offdate,'lat', lat, 'lon', lon, 'elev',elev); %fill in siteStruct with the resulting, sorted information
