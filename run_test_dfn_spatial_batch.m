%RUN_TEST_DFN_SPATIAL_BATCH Entry point for SLURM MATLAB -batch.
try
    out = test_dfn_degradation_spatial_checkpoints(); %#ok<NASGU>
catch ME
    fprintf(2, 'DFN batch run failed:\n%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    rethrow(ME);
end
