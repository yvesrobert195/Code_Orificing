%% User Input
fprintf('*********************************************************************\n')
fprintf('ORIFICE OPTIMIZATION SCRIPT WITH BLANKET TEMPERATURE\n')
fprintf(datestr(now))
fprintf('\n*********************************************************************\n')
User_Input
% INPUT SPECIFICALLY FOR BLANKET
Pb.Constraints.xi_blanket= 40;
Pb.Constraints.xi_power = 40;


%% Initialization
fprintf('*********************************************************************\n')
fprintf('Initialization\n')
fprintf('*********************************************************************\n')
% Assign geometry
fprintf('\tReading geometry\n')
Geometry=geometry(Core);

if strcmpi(Core,'A')
    Geometry.nrings=6; % frontier between power assemblies and blanket
elseif strcmpi(Core,'bnb')
    Geometry.nrings=9;
end

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
Coolant=coolant_properties('sodium',Input,Pb,Geometry);

% Find which assemblies have adjacents and how many adjacent pairs
Input.adjacentAssemblies = findAdjacentAssemblies(Pb.Var.nass, Input.map, Input.lengthQ_original);
Input.nadj = nnz(Input.adjacentAssemblies); % number of adjacent assembly pairs

% Read and create rings
Geometry.rings=find_rings(Input.adjacentAssemblies);

% Create constraints
fprintf('*********************************************************************\n')
fprintf('Constraint matrices\n')
fprintf('*********************************************************************\n')
Pb=make_constraints_blanket(Pb,Input,Coolant,Geometry);

% Solve
fprintf('*********************************************************************\n')
fprintf('Solving with CPLEX\n');
fprintf('*********************************************************************\n')
[Solution.solutionvector, Solution.objval, Solution.status, Solution.output] = cplexmilp(Pb.CPLEX.c, Pb.CPLEX.Aineq, Pb.CPLEX.bineq, Pb.CPLEX.Aeq, Pb.CPLEX.beq, [], [], [], [], [], Pb.CPLEX.ctype,[],Pb.CPLEX.opts);
fprintf('exit status = % i\n', Solution.status);
fprintf('solution time = % f\n', Solution.output.time);

% Post processing and check
if Solution.status==1
    fprintf('*********************************************************************\n')
    fprintf('SOLUTION FOUND : Post-processing and checking errors\n');
    fprintf('*********************************************************************\n')
    Output=post_process_blanket(Solution,Pb,Input,Coolant,Geometry);
else
    fprintf('\n**********************    NO SOLUTION FOUND    **********************\n\n')
end
fprintf('\tSaving workspace\n');
save Last_Problem

fprintf('*********************************************************************\n')
fprintf('END OF SCRIPT\n');
fprintf('*********************************************************************\n')
diary off

try
    movefile clone* ./Clones/;
end

