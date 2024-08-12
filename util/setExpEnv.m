global ISPARALLEL;
%ISPARALLEL = true;
 ISPARALLEL = false;
utilFcnPath = genpath( 'util' );
addpath( utilFcnPath );

if ISPARALLEL
    parpool_least_version = '8.2.0'; % R2013b
    if verLessThan( 'matlab', parpool_least_version )
        % checking to see if my pool is already open
        if matlabpool('size') == 0 
            matlabpool%('open','local',nbCore);
        end
    else
        % To check if a pool is already started
        if isempty(gcp('nocreate'))
           parpool;
        end
    end
end
