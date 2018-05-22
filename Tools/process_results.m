clear

list={'sol_05-21-18_12-27-30_A4var_ref.mat'};%,...
    %'sol_05-18-18_11-55-19_A2var.mat'};

for l=1:length(list)
    names{l}=list{l}(23:end-4);
    load(['../Solutions/' list{l}]);
    for i=1:length(Results)
        Values(i,1)=Results{i,1};
        if ~isempty(Results{i,3}.objval)
            Objval(i,l)=Results{i,3}.objval;
            Time(i,l)=Results{i,3}.output.time;
            Iter(i,l)=Results{i,3}.output.iterations;
        else
            Objval(i,l)=NaN;
            Time(i,l)=NaN;
            Iter(i,l)=NaN;
        end
    end
    Objval(Objval==0)=NaN;
    Time(Time==0)=NaN;
    Iter(Iter==0)=NaN;
    
    eval([names{l} ' = Results;']);
end

figure
plot(Values,Objval);legend(names);title('Objval');
figure
semilogy(Values,Time);legend(names);title('Time');
figure
semilogy(Values,Iter);legend(names);title('Iterations');