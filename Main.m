clear variables;
addpath(genpath('Functions'));
addpath(genpath('Tools'));

delete Last_run.log
diary('Last_run.log')

%% User Input
fprintf('*********************************************************************\n')
fprintf('ORIFICE OPTIMIZATION SCRIPT\n')
fprintf(datestr(now))
fprintf('\n*********************************************************************\n')

fprintf('Reading input\n')
Core='BNB';
%Input.powerDetectorFiles = {'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/base/A_det0.m'};%,...
% '/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/base/A_det1.m'};
Input.powerDetectorFiles = {'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/bnb/BnB_det0.m'};%,...
% '/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/bnb/BnB_det1.m',...
%'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/bnb/BnB_det2.m',...
%'/Users/yvesrobert/Desktop/Stage/Control_rods/try/Results/bnb/BnB_det3.m'};
Input.assemblyPowerThreshold=0.01E6;
Pb.Var.x=logspace(-2,2,1);

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

%% Initialization
fprintf('*********************************************************************\n')
fprintf('Initialization\n')
fprintf('*********************************************************************\n')
% Assign geometry
fprintf('\tReading geometry\n')
Geometry=geometry(Core);

% Reads det files
fprintf('\tReading power data\n')
Input.Q = readQ(Input.powerDetectorFiles);
Input.lengthQ_original = length(Input.Q);
[Input.Q, Input.map] = formMap(Input.Q, Input.assemblyPowerThreshold);
Input.Q_ave = sum(Input.Q,2)/length(Input.powerDetectorFiles); % divide by number of steps to get average assembly power over cycle

fprintf('\tManipulating data\n');
Pb.Var.nass = length(Input.Q_ave); % number of assemblies in problem
Pb.Var.nsteps = length(Input.powerDetectorFiles);
Pb.Var.npossflows = length(Pb.Var.x); % number of possible flowrates specified as data
Pb.Var.nvars = Pb.Var.nass+Pb.Var.nass*Pb.Var.npossflows+Pb.Var.npossflows; % number of total variables

% Build tables of coolant properties and temperature
fprintf('\tBuilding coolant properties and temperature tables\n')
Coolant=coolant_properties('sodium',Input,Pb,Geometry,'c');

% Find which assemblies have adjacents and how many adjacent pairs
Input.adjacentAssemblies = findAdjacentAssemblies(Pb.Var.nass, Input.map, Input.lengthQ_original);
Input.nadj = nnz(Input.adjacentAssemblies); % number of adjacent assembly pairs

% Create constraints
fprintf('*********************************************************************\n')
fprintf('Constraint matrices\n')
fprintf('*********************************************************************\n')
Pb=make_constraints(Pb,Input,Coolant,Geometry);

% Solve
fprintf('*********************************************************************\n')
fprintf('Solving with CPLEX\n');
fprintf('*********************************************************************\n')
[Solution.solutionvector, Solution.objval, Solution.status, Solution.output] = cplexmilp(Pb.CPLEX.c, Pb.CPLEX.Aineq, Pb.CPLEX.bineq, Pb.CPLEX.Aeq, Pb.CPLEX.beq, [], [], [], [], [], Pb.CPLEX.ctype,[],Pb.CPLEX.opts);
fprintf('exit status = % i\n', Solution.status);
fprintf('solution time = % f\n', Solution.output.time);

% Post processing and check
if ~isempty(Solution.solutionvector)
    fprintf('*********************************************************************\n')
    fprintf('Post-processing and checking errors\n');
    fprintf('*********************************************************************\n')
    Output=post_process(Solution,Pb,Input,Coolant,Geometry);
else
    fprintf('\n**********************    NO SOLUTION FOUND    **********************\n\n')
end
fprintf('\tSaving workspace\n');
save Last_Problem

rmpath(genpath('Functions'))
fprintf('*********************************************************************\n')
fprintf('END OF SCRIPT\n');
fprintf('*********************************************************************\n')
diary off