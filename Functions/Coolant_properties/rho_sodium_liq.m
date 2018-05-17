function rho=rho_sodium_liq(T)
A12=1.00423e3;
A13=-0.21390;
A14=1.1046e-5;
rho=A12+A13*T+A14*T.^2;
% rho=850;