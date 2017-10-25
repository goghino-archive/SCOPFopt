% Author: Samuel A. Cruz Alegría.
% Version: 2017-10-13 (YYYY-MM-DD).

% This experiment performs benchmarks for the following cases:
%  1. Number of iterations w.r.t. tolerance value.
function [SUCCESS, INFO] = myExperiment2(mpc, mpopt, savePlots, pathToSavePlots, ...
                                 saveWorkspace, loadWorkspace, ...
                                 workspaceName, untilConvergence)

%% Do computations if necessary.
if loadWorkspace
   disp('Loading workspace');
   load('-mat', workspaceName);
else
   % Fixed contingency set. 
%    contingencies = 10:5:60;
   contingencies = ones(1, 50);
   
   maxTolPowIndex = 2;
   minTolPowIndex = 9;
   step = 1;
   
   % Number of times to repeat the experiment.
   experimentRepetitions = 10;
   
   % All iterations.
   allIters = zeros(experimentRepetitions, length(maxTolPowIndex:step:minTolPowIndex));
   % Tolerance values array.
   tolValuesArray = maxTolPowIndex:step:minTolPowIndex;
   
   if ~untilConvergence
      mpopt.ipopt.opts.max_iter = 5;
   end
   
   for i = 1:experimentRepetitions
      iterWrtTol = [];
      
      % Run SCOPF solver on varying tolerance values.
      % TODO: Set SUCCESS to 0 if any exceptions arise.
      for powIndex = maxTolPowIndex:step:minTolPowIndex
         mpopt.ipopt.opts.tol = 1 * 10^(-powIndex);
         
         [RESULTS, SUCCESS, info] = runscopf(mpc, contingencies, mpopt);
         
         iterWrtTol = [iterWrtTol info.numIter];
      end
      
      allIters(i, :) = iterWrtTol;
   end
end

%% Compute averages.

% Average arrays of collected data.
avgIterWrtTol = zeros(1, size(allIters,2));

for i = 1:size(allIters,2)
   avgIterWrtTol(i) = mean(allIters(:, i));
end

INFO.tolValuesArray = tolValuesArray;
INFO.avgIterWrtTol = avgIterWrtTol;

SUCCESS = 1;

if saveWorkspace
   save(workspaceName);
end

end
