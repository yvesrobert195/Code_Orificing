function f=novendstern(m,T,Pw,pitch,D,H)
Re=Reynolds(Pw,m,T);
M=(1.034/((pitch/D).^0.124)+(29.7.*((pitch/D).^6.94).*Re.^0.086)./((H/D).^2.239)).^0.885;
f_smooth=0.316./(Re.^0.25);
f=M*f_smooth;