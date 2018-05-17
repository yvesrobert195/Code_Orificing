function f=plot_core(corename,nrings,assemblyMap)
%% First part taken from convertSerpentLatticeForHexes.m (maybe can be shortened)
%clearvars -except T_out_hexmap m_hexmap;
if ~exist('assemblyMap')
    assemblyMap = []; %leave blank if you want to plot power instead
end
% assemblyMap = [m_hexmap]; %leave blank if you want to plot power instead

if length(assemblyMap) >= 1 %if assembly map is provided, just plot assembly map
    Q_tab = assemblyMap; %confusing, but simple
    nrows = length(Q_tab);
else %otherwise plot the powers read from corename
    load(corename)
    
    %axially integrate assembly powers
    DETPower=who('DETAssemblyPowerAxial*');
    Q = 0;
    for i=1:length(DETPower)
        Pow=eval(DETPower{i});
        Q =  Q + Pow(:,11);
    end
    totalPower = sum(Q);
    %fission power in each assembly, W
    %        Q_gamma = zeros(size(Q));
    if ~isempty(who('DETAssemblyGammaPowerAxial*'))
        DETGamma=who('DETAssemblyGammaPowerAxial*');
        Q_gamma=0;
        for i=1:length(DETGamma)
            Pow_gamma=eval(DETGamma{i});
            Q_gamma =  Q_gamma + Pow_gamma(:,11);
        end
        %add gamma power into total power
        totalGammaPower = sum(Q_gamma);
        fissionFrac = totalPower/(totalPower+totalGammaPower);
        Q = fissionFrac*Q; %scale total power by the fraction of power from fission, since this is the tally
        Q = Q + Q_gamma; %add in the gamma power
    end
    %     convert Q list to square table
    nrows = sqrt(length(Q));
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
    if ~exist('assemblyMap') || isempty(assemblyMap)
        Hexes_to_Plot(i,3)=Q_tab(rc(i,1),rc(i,2))/10^6;
    else
        Hexes_to_Plot(i,3)=Q_tab(rc(i,1),rc(i,2));
    end
    i = i + 1;
end

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

%% Plotting part
f=figure;
axis equal
colormap(jet(1000))
c=ones(size(colormap));
c(0.25*length(c):0.65*length(c),:)=0;
colorbar
for i=1:length(x)
    col=round((length(c)-1)*(Hexes_to_Plot(i,3)-nanmin(Hexes_to_Plot(:,3)))/(nanmax(Hexes_to_Plot(:,3))-nanmin(Hexes_to_Plot(:,3))));
    if col~=0 && ~isnan(Hexes_to_Plot(i,3))
        p=patch(xhex+x(i),yhex+y(i),Hexes_to_Plot(i,3),'EdgeColor','none'); % make a hexagon at [2i,2j]        ndim=[(dim(1)-AX(1))/Xrange+AX(1)/Xrange (dim(2)-AX(3))/Yrange+AX(3)/Yrange (dim(3)-AX(1))/Xrange (dim(4)-AX(3))/Yrange];
        tx=text(x(i), y(i), num2str(round(Hexes_to_Plot(i,3),2)),'HorizontalAlignment','center','VerticalAlignment','middle','Color',c(col,:), 'FontSize', 7);
        hold on
    end
end
try
    filename=strfind(corename,'/');
    ax = gca;
    ax.Visible = 'off';
    title(corename(filename(end)+1:end),'Interpreter','none','Visible','on')
end