% Author: Samuel A. Cruz Alegría.
% Version: 2017-10-13 (YYYY-MM-DD)
clc;
clear;
close all;

define_constants;

setenv('OMP_NUM_THREADS', '1')

% Can only use one at a time.
% REVIEW: Change this so that you can run experiments for both at once.
usingIPOPT = 0;
usingOptizelle = 1;

% mpopt = mpoption('opf.ac.solver', 'MIPS', 'verbose', 2, 'mips.max_it', 100);

if usingIPOPT
   % Path to save plots in.
   pathToSavePlots = '/Users/samuelcruz/Documents/GitHub/bachelor-project/plots/ipopt/';
   mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);
elseif usingOptizelle
   pathToSavePlots = '/Users/samuelcruz/Documents/GitHub/bachelor-project/plots/optizelle/';
   mpopt = mpoption('opf.ac.solver', 'OPTIZELLE');
end

% load MATPOWER case struct, see help caseformat
mpc = loadcase('case9');
% mpc = loadcase('case118');

runExperiment1 = 1;
runExperiment2 = 0;

%% Experiment 1.
if runExperiment1
   %% Settings for experiment 1.
   % Toggle on/off to save plots.
   savePlots = 0;
   % Toggle on/off to save workspace.
   saveWorkspace = 0;
   % Toggle on/off to load workspace.
   loadWorkspace = 0;
   % Toggle on/off to run until convergence
   untilConvergence = 1;
   
   % Path to save workspace in.
   if untilConvergence
      if usingIPOPT
         workspaceName = ['/Users/samuelcruz/Documents/GitHub/'...
            'bachelor-project/workspace/ipopt_exp1_conv_data.mat'];
      elseif usingOptizelle
         workspaceName = ['/Users/samuelcruz/Documents/GitHub/'...
            'bachelor-project/workspace/optizelle_exp1_conv_data.mat'];
      end
   else
      if usingIPOPT
         workspaceName = ['/Users/samuelcruz/Documents/GitHub/'...
            'bachelor-project/workspace/ipopt_exp1_maxiter_data.mat'];
      elseif usingOptizelle
         workspaceName = ['/Users/samuelcruz/Documents/GitHub/'...
            'bachelor-project/workspace/optizelle_exp1_maxiter_data.mat'];
      end
   end
   
   %% Run experiment 1.
   [SUCCESS, INFO] = myExperiment1(mpc, mpopt, savePlots, pathToSavePlots,...
      saveWorkspace, loadWorkspace, workspaceName, untilConvergence);
   
   %% Experiment 1 plots.
   numContingenciesArray = INFO.numContingenciesArray;
   avgIterWrtNumContingencies = INFO.avgIterWrtNumContingencies;
   avgTimeWrtNumContingencies = INFO.avgTimeWrtNumContingencies;
   stdAllTimes = INFO.stdAllTimes;
   stdAllIters = INFO.stdAllIters;
   
   % Plot average number of iterations vs number of contingencie, with
   % standard deviation.
%    plot(numContingenciesArray, avgIterWrtNumContingencies, '-o');
   errorbar(numContingenciesArray, avgTimeWrtNumContingencies, stdAllIters, '-o');
   title('Average no. of iterations vs No. of contingencies');
   xlabel('No. of contigencies');
   ylabel('Average no. of iterations');
   if savePlots
      if usingIPOPT
         print(strcat(pathToSavePlots, 'ipopt_avgIterAndStdVsNumCont'), '-depsc');
      elseif usingOptizelle
         print(strcat(pathToSavePlots, 'optizelle_avgIterAndStdVsNumCont'), '-depsc');
      end
   end
     
   % Plot average time taken vs number of contingencies, with standard
   % deviation.
   errorbar(numContingenciesArray, avgTimeWrtNumContingencies, stdAllTimes, '-o');
   title('Average time taken vs No. of contingencies');
   xlabel('No. of contigencies');
   ylabel('Average time taken (seconds)');
   if savePlots
      if usingIPOPT
         print(strcat(pathToSavePlots, 'ipopt_avgTimeAndStdInSecsVsNumCont'), '-depsc');
      elseif usingOptizelle
         print(strcat(pathToSavePlots, 'optizelle_avgTimeAndStdInSecsVsNumCont'), '-depsc');
      end
   end
   
   % Plot average time normalized by average number of iterations.
   avgTimePerIter = avgTimeWrtNumContingencies ./ avgIterWrtNumContingencies;
   plot(numContingenciesArray, avgTimePerIter, '-o')
   title('Average time per iteration vs No. of contingencies');
   xlabel('No. of contigencies');
   ylabel('Average time per iteration (seconds)');
   if savePlots
      if usingIPOPT
         print(strcat(pathToSavePlots, 'ipopt_avgTimePerIterInSecsVsNumCont'), '-depsc');
      elseif usingOptizelle
         print(strcat(pathToSavePlots, 'optizelle_avgTimePerIterInSecsVsNumCont'), '-depsc');
      end
   end
end

%% Experiment 2.
if runExperiment2   
   %% Settings for experiment 2.
   % Toggle on/off to save plots.
   savePlots = 0;
   % Toggle on/off to save workspace.
   saveWorkspace = 0;
   % Toggle on/off to load workspace.
   loadWorkspace = 0;
   % Toggle on/off to run until convergence
   untilConvergence = 1;

   % Path to save workspace in.
   if untilConvergence
      if usingIPOPT
         workspaceName = ['/Users/samuelcruz/Documents/GitHub/'...
         'bachelor-project/workspace/ipopt_exp2_conv_data.mat'];
      elseif usingOptizelle
         workspaceName = ['/Users/samuelcruz/Documents/GitHub/'...
         'bachelor-project/workspace/optizelle_exp2_conv_data.mat'];
      end
   else
      if usingIPOPT
         workspaceName = ['/Users/samuelcruz/Documents/GitHub/'...
            'bachelor-project/workspace/ipopt_exp2_maxiter_data.mat'];
      elseif usingOptizelle
         workspaceName = ['/Users/samuelcruz/Documents/GitHub/'...
            'bachelor-project/workspace/optizelle_exp2_maxiter_data.mat'];
      end
   end
   
   %% Run experiment 2.
   [SUCCESS, INFO] = myExperiment2(mpc, mpopt, savePlots, pathToSavePlots,...
      saveWorkspace, loadWorkspace, workspaceName, untilConvergence);
   
   %% Experiment 2 plots
   % Plot average number of iterations vs number of contingencies.
   plot(INFO.tolValuesArray, INFO.avgIterWrtTol, '-o');
   title('Average no. of iterations vs Tolerance value');
   xlabel('-log(tolerance value)');
   ylabel('Average no. of iterations');
   if savePlots
      if usingIPOPT
         print(strcat(pathToSavePlots, 'ipopt_avgIterVsTol'), '-depsc');
      elseif usingOptizelle
         print(strcat(pathToSavePlots, 'optizelle_avgIterVsTol'), '-depsc');
      end
   end
end