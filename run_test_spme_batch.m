function run_test_spme_batch()
%RUN_TEST_SPME_BATCH Cluster-friendly wrapper for long SPMe test.

    project_root = fileparts(mfilename('fullpath'));
    addpath(project_root);
    addpath(fullfile(project_root, 'tests'));

    results_dir = fullfile(project_root, 'results');
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    ts = datestr(now, 'yyyymmdd_HHMMSS');
    diary_file = fullfile(results_dir, ['matlab_diary_' ts '.txt']);
    diary(diary_file);
    diary on;

    logmsg('Job started');
    logmsg(['MATLAB diary: ' diary_file]);

    try
        logmsg('Running test_spme_ramadass_vs_yang...');
        out = test_spme_ramadass_vs_yang();

        final_file = fullfile(results_dir, ['test_spme_ramadass_vs_yang_' ts '.mat']);
        logmsg(['Saving final MAT file: ' final_file]);
        save(final_file, 'out', '-v7.3');

        logmsg('Job finished successfully');
    catch ME
        fprintf(2, '[%s] ERROR: %s\n', datestr(now, 31), ME.message);
        for k = 1:numel(ME.stack)
            fprintf(2, '  In %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
        end

        err_file = fullfile(results_dir, ['test_spme_ramadass_vs_yang_ERROR_' ts '.mat']);
        save(err_file, 'ME', '-v7.3');
        logmsg(['Saved error snapshot: ' err_file]);

        diary off;
        rethrow(ME);
    end

    diary off;
end

function logmsg(message)
    fprintf('[%s] %s\n', datestr(now, 31), message);
end
