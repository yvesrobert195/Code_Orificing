function old_plot_vect(corename,nass,vect)
% %% take a core
% corename='./twor/Core_0.mat';
% nass=523;
%% read a det file
load(corename)
Q = [];
DETPower=who('DETAssemblyPowerAxial*');
Q_step = 0;
for i=1:length(DETPower)
    Pow=eval(DETPower{i});
    Q_step =  Q_step + Pow(:,11);
end

DETGamma=who('DETAssemblyGammaPowerAxial*');
Q_gamma_step=0;
for i=1:length(DETGamma)
    Pow_gamma=eval(DETGamma{i});
    Q_gamma_step =  Q_gamma_step + Pow_gamma(:,11);
end
%     Q_step = DETAssemblyPowerAxial1(:,11)+DETAssemblyPowerAxial2(:,11)+DETAssemblyPowerAxial3(:,11)+DETAssemblyPowerAxial4(:,11)+DETAssemblyPowerAxial5(:,11)+DETAssemblyPowerAxial6(:,11)+DETAssemblyPowerAxial7(:,11)+DETAssemblyPowerAxial8(:,11)+DETAssemblyPowerAxial9(:,11)+DETAssemblyPowerAxial10(:,11)+DETAssemblyPowerAxial11(:,11)+DETAssemblyPowerAxial12(:,11)+DETAssemblyPowerAxial13(:,11)+DETAssemblyPowerAxial14(:,11)+DETAssemblyPowerAxial15(:,11)+DETAssemblyPowerAxial16(:,11)+DETAssemblyPowerAxial17(:,11)+DETAssemblyPowerAxial18(:,11)+DETAssemblyPowerAxial19(:,11)+DETAssemblyPowerAxial20(:,11)+DETAssemblyPowerAxial21(:,11)+DETAssemblyPowerAxial22(:,11)+DETAssemblyPowerAxial23(:,11)+DETAssemblyPowerAxial24(:,11); %fission power in each assembly, W
%     Q_gamma_step = DETAssemblyGammaPowerAxial1(:,11)+DETAssemblyGammaPowerAxial2(:,11)+DETAssemblyGammaPowerAxial3(:,11)+DETAssemblyGammaPowerAxial4(:,11)+DETAssemblyGammaPowerAxial5(:,11)+DETAssemblyGammaPowerAxial6(:,11)+DETAssemblyGammaPowerAxial7(:,11)+DETAssemblyGammaPowerAxial8(:,11)+DETAssemblyGammaPowerAxial9(:,11)+DETAssemblyGammaPowerAxial10(:,11)+DETAssemblyGammaPowerAxial11(:,11)+DETAssemblyGammaPowerAxial12(:,11)+DETAssemblyGammaPowerAxial13(:,11)+DETAssemblyGammaPowerAxial14(:,11)+DETAssemblyGammaPowerAxial15(:,11)+DETAssemblyGammaPowerAxial16(:,11)+DETAssemblyGammaPowerAxial17(:,11)+DETAssemblyGammaPowerAxial18(:,11)+DETAssemblyGammaPowerAxial19(:,11)+DETAssemblyGammaPowerAxial20(:,11)+DETAssemblyGammaPowerAxial21(:,11)+DETAssemblyGammaPowerAxial22(:,11)+DETAssemblyGammaPowerAxial23(:,11)+DETAssemblyGammaPowerAxial24(:,11); %gamma power in each assembly, W
%
%add gamma power into total power
totalPower = sum(Q_step);
totalGammaPower = sum(Q_gamma_step);
fissionFrac = totalPower/(totalPower+totalGammaPower);
Q_step = fissionFrac*Q_step; %scale total power by the fraction of power from fission, since this is the tally
Q_step = Q_step + Q_gamma_step; %add in the gamma power

%add assembly powers in current depletion step to running total
Q(:,end+1) = Q_step;
%% convert it to a map
 nrows = sqrt(length(Q));
 nrings=nrows/2+1;
    j = 1;
    while j < nrows+1
        i = 1;
        while i < nrows+1
            Q_tab(j,i) = Q(1);
            Q(1) = [];
            i = i + 1;
        end
        j = j + 1;
    end

if (nrings-1)*2 > nrows
    display('***error!')
    display('***you have entered more rings than there are in the power matrix')
end

rc = []; %row, column
rp = []; %ring, position

i = 1; %assembly count
%special case of center assembly
rc(1,1) = (nrows+1)/2; rc(1,2) = (nrows+1)/2;
rp(1,1) = 1; rp(1,2) = 1;
i = i + 1;

%all rings other than center assembly

for ring = 2:nrings
    rowsignal = []; %used to store the pattern for each row of the matrix, which looks like a negative sin wave
    j = 0;
    while j < ring
        rowsignal(end+1) = -j;
        j = j + 1;
    end
    j = 1;
    while j < ring
        rowsignal(end+1) = -(ring-1);
        j = j + 1;
    end
    j = -(ring-2);
    while j < ring
        rowsignal(end+1) = j;
        j = j + 1;
    end
    j = 1;
    while j < ring
        rowsignal(end+1) = ring-1;
        j = j + 1;
    end
    j = ring-2;
    while j > 0
        rowsignal(end+1) = j;
        j = j - 1;
    end
    
    columnsignal = [rowsignal(2*(ring-1)+1:end), rowsignal(1:2*(ring-1))]; %the column signal is simply a shifted version of the row signal
    
    rc = [rc; rowsignal'+(nrows+1)/2,columnsignal'+(nrows+1)/2]; %populate the row/column array
    
    pos = 1; %position in ring count
    while pos < 6*(ring-1)+1 %populate the ring/position array
        rp(i,1) = ring;
        rp(i,2) = pos;
        i = i + 1;
        pos = pos + 1;
    end
end

Hexes_to_Plot=nan*ones(length(rp),3);
%'Ring, Position, Color, Value\n');
i = 1;
while i < length(rp)+1
    Hexes_to_Plot(i,1) =rp(i,1);
    Hexes_to_Plot(i,2) =rp(i,2);
    H=Q_tab(rc(i,1),rc(i,2));
    Hexes_to_Plot(i,3)=H;
    i = i + 1;
end
h=Hexes_to_Plot(:,3);
threshold=min(h);
while sum(~isnan(h))>nass
    h(h<threshold)=NaN;
    threshold=threshold+1;
end
Hexes_to_Plot(:,3)=h;
%% read vector (temp, flowrate, ...)
n=1;
for i=1:size(Hexes_to_Plot(:,3))
    if ~isnan(Hexes_to_Plot(i,3))
        MAP(i,1)=vect(n);
        n=n+1;
    else
        MAP(i,1)=NaN;
    end
end        
Hexes_to_Plot(:,4)=MAP;

%% plot results
R=1;
r=R*sqrt(3)/2;
xdir=r*2*cos(pi/3.*(0:5));
ydir=r*2*sin(pi/3.*(0:5));
xdir=[xdir(3:end) xdir(1:2)];
ydir=[ydir(3:end) ydir(1:2)];
%nw w sw se e ne

xhex=R*cos(pi/6+pi/3.*(0:5)); % x-coordinates of the vertices
yhex=R*sin(pi/6+pi/3.*(0:5)); % y-coordinates of the vertices

x(1)=-r;
y(1)=3/2*R;
for i=2:nrings
    x(end+1)=x(1)+(i-1)*r*2;
    y(end+1)=y(1);
    for j=2:(i-1)*6
        x(end+1)=x(end)+xdir(ceil((j-1)/(i-1)));
        y(end+1)=y(end)+ydir(ceil((j-1)/(i-1)));
    end
    i = i + 1;
end

% Plotting part
f=figure;
ax=axis;
axis equal
colormap(jet(1000))
c=ones(size(colormap));
c(0.25*length(c):0.65*length(c),:)=0;
colorbar
for i=1:length(x)
    col=round((length(c)-1)*(Hexes_to_Plot(i,4)-nanmin(Hexes_to_Plot(:,4)))/(nanmax(Hexes_to_Plot(:,4))-nanmin(Hexes_to_Plot(:,4))));
    if col==0
        col=1;
    end
    if col~=0 && ~isnan(Hexes_to_Plot(i,4))
      p=patch(xhex+x(i),yhex+y(i),Hexes_to_Plot(i,4),'EdgeColor','none'); % make a hexagon at [2i,2j]
        tx=text(x(i), y(i), num2str(round(Hexes_to_Plot(i,4),2)),'HorizontalAlignment','center','VerticalAlignment','middle','Color',c(col,:),'FontUnits', 'Normalized','FontSize',2e-2);
        hold on
    end
end
try
filename=strfind(corename,'/');
ax = gca;
ax.Visible = 'off';
title(corename(filename(end)+1:end),'Interpreter','none','Visible','on')
end