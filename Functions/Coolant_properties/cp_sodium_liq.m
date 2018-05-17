function cp=cp_sodium_liq(T)
%Cl in J/kg-K
A28=7.3898e5;
A29=3.154e5;
A30=1.1340e3;
A31=-2.2153e-1;
A32=1.1156e-4;
Tc=2503.3;
cp=A28./((Tc-T).^2)+A29./(Tc-T)+A30+A31.*(Tc-T)+A32.*(Tc-T).^2;
% cp=1272;
end