function run_tests()

disp('-> run_tests.m')
warning  on

ANTELOPE = getenv('ANTELOPE');
if isempty(ANTELOPE)
    ANTELOPE='/opt/antelope/5.7';
end

% Add GISMO
%NEWGISMOPATH = '/raid/apps/src/GISMO';
NEWGISMOPATH = '/home/a/akfarrell/Uturuncu/GISMO';
addpath(NEWGISMOPATH)
startup_GISMO
disp(sprintf('GISMO toolbox added OK from %s',NEWGISMOPATH));
disp('  - Find out more by visiting http://github.com/geoscience-community-codes/GISMO/wiki)')


% Add Antelope
pathToSetupFile = fullfile(ANTELOPE,'setup.m');
if exist(pathToSetupFile)
	disp(sprintf('- running %s',pathToSetupFile'));
	try
        	run(pathToSetupFile);
		disp('Antelope setup.m OK');
		disp('(Find out more with help antelope or help antelope/examples)')
	catch
		addpath(ANTELOPE);
		setup
	end
else
	warning(sprintf('- %s not found',pathToSetupFile'));
        warning('Antelope toolbox not added');
end


try 
    fprintf('Antelope toolbox for MATLAB is ');
    if ~admin.antelope_exists()
        fprintf('NOT')
    end
    fprintf(' available\n');
    
catch
    warning('admin.antelope_exists() crashed. Antelope NOT available')
end

if exist('/raid/apps/antelope')
    try
    dbptr = dbopen('/raid/apps/antelope/data/db/demo/demo', 'r');
    if isnumeric(dbptr.database) & dbptr.database >= 0
	disp('demo database opens - Antelope Toolbox for MATLAB works')
    else
	warning('demo database does not open - Antelope Toolbox for MATLAB not working')
    end
    catch
        warning('Problem with Antelope Toolbox for MATLAB MEX-files')
    end
else
    warning('Directory /raid/apps/antelope NOT available')
end


disp('<- run_tests.m')
