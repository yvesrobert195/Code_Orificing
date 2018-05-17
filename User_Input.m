fprintf('Reading input\n')
Core='A';
Input.powerDetectorFiles = {'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/base/A_det0.m',...
'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/base/A_det1.m'};
%Input.powerDetectorFiles = {'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/bnb/BnB_det0.m'};%,...
% '/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/bnb/BnB_det1.m',...
%'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/bnb/BnB_det2.m',...
%'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/bnb/BnB_det3.m'};
Input.assemblyPowerThreshold=0.01E6;
Pb.Var.x=logspace(-2,2,45);

% User constraints
Pb.Constraints.T_out_bar = 510; % perfectly mixed coolant outlet plenum temperature, C
Pb.Constraints.T_out_bar_tol=5; % tolerance on mixed outlet temperature
Pb.Constraints.v_max = 12.0; % maximum coolant velocity allowed in assembly, m/s, limit taken from Qvist et al
Pb.Constraints.xi = 40.0; % maximum outlet temperature difference between adjacent assemblies, C
Pb.Constraints.dT_max = 210;
Pb.Constraints.dP_max = 1E6; % maximum allowable pressure drop over core, Pa, limit taken from Qvist et al

% CPLEX options
Pb.CPLEX.opts=cplexoptimset('display','on'); % Option to display iterations ('iter','on','off')
Pb.CPLEX.opts.exportmodel = 'model.lp'; % Name of the saved model
