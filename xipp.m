function O=xipp(S,P,I,C,G)
nass=P.Var.nass;
npossflows=P.Var.npossflows;
nsteps=P.Var.nsteps;
nvars=P.Var.nvars;

Q=I.Q;

cp=C.heat_capacity;
Omega=C.T_gradient;
T_in=C.T_inlet;

% Calculate quantities
fprintf('\tCreating Solution\n')
O.m = S.solutionvector(1:P.Var.nass);
for i=1:nass
    for j=1:npossflows
        O.delta(i,j)=S.solutionvector(nass+(i-1)*npossflows+j);
    end
end
O.xi=S.solutionvector(end);

for k=1:nsteps
    for i=1:nass
        O.Tout(i,k) = T_in + sum(O.delta(i,:).*Omega(i,:,k));
        O.density(i,k) = sum(C.density(i,:,k).*O.delta(i,:));
        O.heat_capacity(i,k) = sum(cp(i,:,k).*O.delta(i,:));
        O.P_gradient(i,k) = sum(C.P_gradient(i,:,k).*O.delta(i,:));
        alpha(i,k)=Q(i,k)/sum(O.delta(i,:).*cp(i,:,k));
    end
end

adjT=NaN*ones(nass,6,nsteps);
for k=1:nsteps
    for i=1:nass
        for j=1:6
            if I.adjacentAssemblies(i,j)~=0
                O.adjT(i,j,k)=abs(O.Tout(i,k)-O.Tout(I.adjacentAssemblies(i,j),k));
            end
        end
    end
end

%check max outlet temp, adj temp and mixed outlet temp constraints
fprintf('\tChecking errors\n')
sum3=0;
er=0;
for k = 1:nsteps
    % Constraint 1
    if sum(O.Tout(:,k) > C.T_inlet+P.Constraints.dT_max) > 0
        fprintf('\terror: an assembly in step %i violates the outlet temperature constraint\n', k);
        er=er+1;
    end
    % Constraint 2
    if sum((C.T_inlet + sum(alpha(G.rings<=P.Constraints.rings_outlet,k))/sum(O.m(G.rings<=P.Constraints.rings_outlet))) > P.Constraints.T_out_bar+P.Constraints.T_out_bar_tol) > 0
        fprintf('\terror: mixed outlet temperature violates constraint in step %i\n', k);
            er=er+1;
    elseif sum((C.T_inlet + sum(alpha(G.rings<=P.Constraints.rings_outlet,k))/sum(O.m(G.rings<=P.Constraints.rings_outlet))) < P.Constraints.T_out_bar-P.Constraints.T_out_bar_tol) > 0
        fprintf('\terror: mixed outlet temperature violates constraint in step %i\n', k);
            er=er+1;
    end
    % ini Constraint 3
    for i=1:nass
        sum3=sum3+(O.m(i)/(O.density(i,k)*G.Assembly.flow_area)>P.Constraints.v_max);
    end
    % Constraint 4
    if max(max(O.adjT(G.rings<G.nrings,:,k)))>O.xi
        fprintf('\terror: Two adjacent power assemblies violate the maximum temperature gradient in step %i\n', k);
            er=er+1;
    end
    % Constraint 9
    if max(O.P_gradient(:,k))>P.Constraints.dP_max
        fprintf('\terror: pressure loss in one assembly violates pressure gradient constraint in step %i\n', k);
            er=er+1;
    end
end

% Constraint 3
if sum3 > 0
    fprintf('\terror: an assembly violates the maximum velocity constraint\n');
        er=er+1;
end

fprintf('\tChecked : ')
if er==0
  fprintf('No error found\n');
else
  fprintf('\tWARNING : %i errors on constraints found during post-processing\n', er);
end

  