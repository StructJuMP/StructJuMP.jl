% Saves the nice output of Matpower (to have it as a comparison for
% StructJuMP)

cases = {'case9', 'case30', 'case118', ...
    'case57', 'case300', 'case1354pegase.m', 'case2869pegase.m', ...
    'case2736sp', 'case2737sop', 'case2383wp', 'case2746wop', 'case2746wop', ...
    'case3012wp', 'case3120sp', 'case3375wp', 'case9241pegase'};

for i=1:length(cases)
   fprintf('running %s ...', cases{i});
   runopf(cases{i}, mpoption(), [cases{i} '.result']) 
   fprintf(' done!\n');
end