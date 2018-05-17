function c=coolant_properties(coolant,I,P,G,constant)
x=P.Var.x;
Q=I.Q;
A_flow=G.Assembly.flow_area;
K=G.Assembly.loss_coefficient;
D_h= G.Assembly.hydraulic_diameter;
H=G.Rod.height;
L=H;

c.density=NaN*ones(size(Q,1),length(x),size(Q,2));
c.P_gradient=NaN*ones(size(Q,1),length(x),size(Q,2));
c.heat_capacity=NaN*ones(size(Q,1),length(x),size(Q,2));
c.T_outletC=NaN*ones(size(Q,1),length(x),size(Q,2));
c.T_outletK=NaN*ones(size(Q,1),length(x),size(Q,2));
c.T_gradient=NaN*ones(size(Q,1),length(x),size(Q,2));

c.T_inlet=355;

if coolant=='sodium'
    
    if exist('constant') && constant=='c'
        c.density(:,:,:)=850;
        c.friction_factor(:,:,:)=0.021132;
        c.heat_capacity(:,:,:)=1272;
        for k=1:size(Q,2)
            c.T_gradient(:,:,k) = Q(:,k)./(x.*c.heat_capacity(:,:,k));
        end
        c.T_outletK(:,:,:)=c.T_gradient(:,:,:)+c.T_inlet+273.15;
        
    else
        f=@(T,q,m) (cp_sodium_liq(T).*(T-(c.T_inlet+273.15))-q/m);
        for k=1:size(Q,2)
            for i=1:size(Q,1)
                for j=1:length(x)
                    fun = @(T) f(T,Q(i,k),x(j));
                    c.T_outletK(i,j,k)=fzero(fun,P.Constraints.T_out_bar);
                    c.heat_capacity(i,j,k)=cp_sodium_liq(c.T_outletK(i,j,k));
                    c.density(i,j,k)=rho_sodium_liq(c.T_outletK(i,j,k));                
                    c.T_gradient(i,j,k) = Q(i,k)./(x(j)*c.heat_capacity(i,j,k));
                end
                
            end
        end
    end
    for k=1:size(Q,2)
        for i=1:size(Q,1)
            for j=1:length(x)
                c.T_outletC(i,j,k)=c.T_outletK(i,j,k)-273.15;
                f=novendstern(x(j),c.T_outletK(i,j,k),G.Assembly.wetted_perimeter,G.Rod.pitch,G.Rod.rod_diameter,G.Rod.height);
                c.P_gradient(i,j,k)=x(j)^2*(K/(2*c.density(i,j,k)*A_flow^2)+f*L/(D_h*2*c.density(i,j,k)*A_flow^2))+c.density(i,j,k)*9.81*H; 
            end
        end
    end
else
    error('Please specify a valid coolant')
end
end