function O=post_process(S,P,I,C,G)
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
for j=1:npossflows
    O.beta(j,1)=S.solutionvector(nass+nass*npossflows+j);
end

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
for k = 1:nsteps
    if sum(O.Tout(:,k) > C.T_inlet+P.Constraints.dT_max) > 0
        fprintf('\terror: an assembly in step %i violates the outlet temperature constraint\n', k);
    end
    if max(max(O.adjT(:,:,k)))>P.Constraints.xi
        fprintf('\terror: Two adjacent assemblies violate the maximum temperature gradient in step %i\n', k);
    end
    if sum((C.T_inlet + sum(alpha(:,k))/sum(O.m)) > P.Constraints.T_out_bar+P.Constraints.T_out_bar_tol) > 0
        fprintf('\terror: mixed outlet temperature violates constraint in step %i\n', k);
    elseif sum((C.T_inlet + sum(alpha(:,k))/sum(O.m)) < P.Constraints.T_out_bar-P.Constraints.T_out_bar_tol) > 0
        fprintf('\terror: mixed outlet temperature violates constraint in step %i\n', k);
    end
end

sum3=0;
for k=1:nsteps
    for i=1:nass
        sum3=sum3+(O.m(i)/(O.density(i,k)*G.Assembly.flow_area)>P.Constraints.v_max);
    end
end
if sum3 > 0
    fprintf('\terror: an assembly violates the maximum velocity constraint\n');
end

fprintf('\tChecked\n')