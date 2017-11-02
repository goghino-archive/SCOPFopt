% Author: Samuel A. Cruz Alegría.
% Version: 2017-10-13 (YYYY-MM-DD)
clc;
clear;
close all;

define_constants;

setenv('OMP_NUM_THREADS', '1')
% mpopt = mpoption('opf.ac.solver', 'MIPS', 'verbose', 2, 'mips.max_it', 100);
mpopt = mpoption('opf.ac.solver', 'Optizelle', 'verbose', 2);
% mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);

% load MATPOWER case struct, see help caseformat
mpc = loadcase('case118');

contingencies = [1];
[RESULTS, SUCCESS, info] = optizellescopf_solver(mpc, contingencies, mpopt);