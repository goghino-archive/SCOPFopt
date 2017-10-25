% Author: Samuel A. Cruz Alegría.
% Version: 2017-10-13 (YYYY-MM-DD).

% This experiment performs benchmarks for the following cases:
%  1. Number of iterations w.r.t. number of contingencies.
%  2. Time taken by the algorithm w.r.t. number of contingencies.
function [SUCCESS, INFO] = myExperiment1(mpc, mpopt, savePlots, pathToSavePlots, ...
                                 saveWorkspace, loadWorkspace, ...
                                 workspaceName, untilConvergence)
                              
SUCCESS = 1;
                              
if loadWorkspace
   disp('Loading workspace');
   load('-mat', workspaceName);
else
   %% Do computations if necessary.
   contingencies = [];
   totalNumContingencies = 100;
   
   % Starting value for number of contingencies.
   lower = 0;
   % Step with which to increase number of contingencies.
   step = 5;
   
   % Number of times to repeat the experiment.
   experimentRepetitions = 10;
   
   % Results of each repetition of the experiment.
   allTimes = zeros(experimentRepetitions, length(lower:step:totalNumContingencies));
   allIters = zeros(experimentRepetitions, length(lower:step:totalNumContingencies));
   
   % Set tolerance value for IPOPT.
   mpopt.ipopt.opts.tol = 1e-6;
   
   if ~untilConvergence
      mpopt.ipopt.opts.max_iter = 5;
   end
   
   for i = 1:experimentRepetitions
      numContingenciesArray = [];
      
      % Number of iterations with respect to number of contingencies.
      iterWrtNumContingencies = [];
      
      % Time taken (in seconds) with respect to number of contingencies.
      timeWrtNumContingencies = [];
      
      % Run SCOPF solver on varying contingency sets.
      % TODO: Set SUCCESS to 0 if any exceptions arise.
      for numContingencies = lower:step:totalNumContingencies
         %         contingencies = lower:(lower + numContingencies - 1);
         contingencies = ones(1, numContingencies);
         
         [RESULTS, SUCCESS, info] = runscopf(mpc, contingencies, mpopt);
         
         numContingenciesArray = [numContingenciesArray numContingencies];
         iterWrtNumContingencies = [iterWrtNumContingencies info.numIter];
         timeWrtNumContingencies = [timeWrtNumContingencies info.overallAlgorithm];
      end
      
      % Add results for current instance of the test.
      allIters(i, :) = iterWrtNumContingencies;
      allTimes(i, :) = timeWrtNumContingencies;
   end
end

%% Compute averages.
% Number of samples.
numSamples = size(allIters, 2);

% Average arrays of collected data.
avgIterWrtNumContingencies = zeros(1, numSamples);
avgTimeWrtNumContingencies = zeros(1, numSamples);

for i = 1:size(allIters,2)
   avgIterWrtNumContingencies(i) = mean(allIters(:, i));
   avgTimeWrtNumContingencies(i) = mean(allTimes(:, i));
end

% Standard deviation in time from each experiment repetition.
stdAllTimes = std(allTimes);

% Standard deviation in number of iterations from each experiment repetition.
stdAllIters = std(allIters);

% Deliver necessary information.
INFO.numContingenciesArray = numContingenciesArray;
INFO.avgIterWrtNumContingencies = avgIterWrtNumContingencies;
INFO.avgTimeWrtNumContingencies = avgTimeWrtNumContingencies;
INFO.stdAllTimes = stdAllTimes;
INFO.stdAllIters = stdAllIters;


% if saveWorkspace
%    disp('Saving workspace...');
%    save(workspaceName);
% end

end
