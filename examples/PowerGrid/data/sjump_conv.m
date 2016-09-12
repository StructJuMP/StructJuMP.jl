% Saves the MATPOWER cases data into tabular data
%
% Usage:
% Run this in the MATPOWER directory containing the cases and specify which
% cases to save in the 'cases' cell array below.
%
% Remarks: 
% The MATPOWER format is a slightly different than the IEEE formats
% (http://www.ee.washington.edu/research/pstca/formats/pti.txt). MATPOWER
% contains fields not included in the PTI data; it also does not use some 
% fields present in PTI (MATPOWER's caseformat.m).

clear all
%cases = {'case9.m', 'case30.m', 'case118.m'};
cases = {'case9.m', 'case30.m', 'case118.m', ...
    'case57.m', 'case300.m', 'case1354pegase.m', 'case2869pegase.m', ...
    'case2736sp.m', 'case2737sop.m', 'case2383wp.m', 'case2746wop.m', 'case2746wop.m', ...
    'case3012wp.m', 'case3120sp.m', 'case3375wp.m', 'case9241pegase.m'}
%cases = {'case1354pegase.m'};
for i=1:length(cases)
    [~,name]=fileparts(cases{i});
    fprintf('converting %s\n', cases{i});
    f=str2func(name);
    mpc=f();
    
    %if size(mpc.bus,2)~=13
    %   msg = sprintf('truncating bus info to only 13 columns for %s', cases{i});
    %   warning(msg)
    %   mpc.bus = mpc.bus(:,1:13);
    %end
    dlmwrite([name, '.bus'], mpc.bus, 'delimiter', '\t', 'precision', '%.8f');
    dlmwrite([name, '.branch'], mpc.branch, 'delimiter', '\t', 'precision', '%.8f');
    dlmwrite([name, '.gen'], mpc.gen, 'delimiter', '\t', 'precision', '%.8f');
    dlmwrite([name, '.gencost'], mpc.gencost, 'delimiter', '\t', 'precision', '%.8f');
end