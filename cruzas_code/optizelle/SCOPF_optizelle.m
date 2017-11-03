% Author: Samuel A. Cruz Alegría.
% Version: 2017-10-13 (YYYY-MM-DD)
clc;
clear;
close all;

setenv('OMP_NUM_THREADS', '1')

% load MATPOWER case struct, see help caseformat
mpc = loadcase('case118');

contingencies = [1, 1];
optizellescopf_solver(contingencies);