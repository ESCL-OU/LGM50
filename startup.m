function startup()
%STARTUP Project startup hook for LGM50-main.
% Run manually from this folder, or add this folder to MATLAB user path.

    casadi_dir = '/Users/yuichikajiura/Applications/casadi-3.7.2-osx_arm64-matlab2018b';
    project_root = fileparts(mfilename('fullpath'));
    tests_dir = fullfile(project_root, 'tests');

    if isfolder(casadi_dir)
        addpath(casadi_dir);
    else
        warning('startup:missingCasadiPath', ...
            'CasADi path not found: %s', casadi_dir);
    end

    if isfolder(tests_dir)
        addpath(tests_dir);
    else
        warning('startup:missingTestsPath', ...
            'Tests path not found: %s', tests_dir);
    end
end
