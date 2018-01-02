% Author: Samuel A. Cruz Alegría.
% Version: 2018-01-01 (YYYY-MM-DD)
clc;
clear;
close all;

define_constants;

setenv('OMP_NUM_THREADS', '1')

% Can only use one at a time.
% TODO: Change this so that you can run experiments for both at once.
usingIPOPT = 0;
usingOptizelle = 1;

if usingIPOPT
   mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);
elseif usingOptizelle
   mpopt = mpoption('opf.ac.solver', 'OPTIZELLE');
end

theCase = 'case9';
fprintf('Testing with %s\n', theCase);

% load MATPOWER case struct, see help caseformat
mpc = loadcase(theCase);

% cont = [2,3,5,6,8,9];
cont = [];

if usingIPOPT
   mpopt.ipopt.opts.tol = 1e-6;
end

for c = 0:length(cont)
   if (c == 0)
      subcont = [];
   else
      subcont = cont(1:c);
   end
   
   [RESULTS, SUCCESS, info] = runscopf(mpc, cont, mpopt);
end