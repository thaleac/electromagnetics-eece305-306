function results = runTests(varargin)
%RUNTESTS  Run the EECE 306 toolbox test suite and print a summary.
%
%   RUNTESTS runs every file matching tests/test_lab*.m, in order, and
%   prints one line per test file plus a summary. Call it from the top of
%   your EECE306-Toolbox folder.
%
%   RUNTESTS('lab', 5) runs only tests/test_lab05.m.
%   RUNTESTS('verbose', true) prints the full error for each failure.
%   R = RUNTESTS(...) also returns a struct array with the results.
%
%   A test file is any script that runs to completion without throwing.
%   Use em.test.assertClose inside it to check values. A test file that
%   throws is reported as a failure, and the message is shown.
%
%   Test scripts execute inside this function's workspace, so the
%   internal state below carries an rt prefix that no reasonable test
%   variable will collide with. Do not modify this file.
%
%   EECE 306, Howard University, Fall 2026.

rt_p__ = inputParser;
rt_p__.addParameter('lab', [], @(x) isnumeric(x));
rt_p__.addParameter('verbose', false, @islogical);
rt_p__.parse(varargin{:});
rt_onlyLab__ = rt_p__.Results.lab;
rt_verbose__ = rt_p__.Results.verbose;

rt_here__ = fileparts(mfilename('fullpath'));
rt_testDir__ = fullfile(rt_here__, 'tests');

if ~isfolder(rt_testDir__)
    error('runTests:noTestFolder', ...
        'No tests folder found at %s. Create it and add test_lab01.m.', rt_testDir__);
end

if isempty(rt_onlyLab__)
    rt_listing__ = dir(fullfile(rt_testDir__, 'test_lab*.m'));
else
    rt_listing__ = dir(fullfile(rt_testDir__, sprintf('test_lab%02d.m', rt_onlyLab__)));
    if isempty(rt_listing__)
        error('runTests:noSuchLab', 'No test file for lab %d.', rt_onlyLab__);
    end
end

if isempty(rt_listing__)
    fprintf('No test files found in %s\n', rt_testDir__);
    results = struct([]);
    return
end

[~, rt_order__] = sort({rt_listing__.name});
rt_listing__ = rt_listing__(rt_order__);
rt_names__ = {rt_listing__.name};

% Make sure the toolbox itself is reachable, without putting +em on the path.
addpath(rt_here__);

fprintf('\n');
fprintf('EECE 306 toolbox test suite\n');
fprintf('%s\n', repmat('-', 1, 62));

rt_nPass__ = 0;
rt_nFail__ = 0;
rt_results__ = struct('name', {}, 'passed', {}, 'seconds', {}, 'message', {});

for rt_k__ = 1:numel(rt_names__)
    rt_stem__ = rt_names__{rt_k__}(1:end-2);
    rt_t0__ = clock();
    rt_passed__ = true;
    rt_msg__ = '';
    try
        run(fullfile(rt_testDir__, rt_names__{rt_k__}));
    catch rt_err__
        rt_passed__ = false;
        rt_msg__ = rt_err__.message;
    end
    try
        rt_dt__ = etime(clock(), rt_t0__);
    catch
        rt_dt__ = NaN;   % a test clobbered the timer state, report NaN
    end

    if rt_passed__
        rt_nPass__ = rt_nPass__ + 1;
        fprintf('  PASS  %-18s  %6.2f s\n', rt_stem__, rt_dt__);
    else
        rt_nFail__ = rt_nFail__ + 1;
        fprintf('  FAIL  %-18s  %6.2f s   %s\n', rt_stem__, rt_dt__, ...
            strtok(rt_msg__, newline));
        if rt_verbose__
            fprintf('        %s\n', rt_msg__);
        end
    end

    rt_results__(end+1) = struct('name', rt_stem__, 'passed', rt_passed__, ...
        'seconds', rt_dt__, 'message', rt_msg__); %#ok<AGROW>
end

fprintf('%s\n', repmat('-', 1, 62));
fprintf('  %d passed, %d failed, %d total\n', rt_nPass__, rt_nFail__, ...
    rt_nPass__ + rt_nFail__);
if rt_nFail__ == 0
    fprintf('  All tests passing.\n');
else
    fprintf('  Re-run with runTests(''verbose'', true) for full messages.\n');
end
fprintf('\n');

if nargout > 0
    results = rt_results__;
end
end
