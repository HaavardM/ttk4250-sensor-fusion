clear all; close all;
%% Parameters
showplt_gnss = true;
showplt_estimates = true;
showplt_state_errors = true;
showplt_error_distance = false;
showplt_nees = true;
showplt_boxplot = false;
%% Run
run_INS_simulated;
run_INS_real;
