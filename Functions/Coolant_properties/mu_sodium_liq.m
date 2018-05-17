function mu=mu_sodium_liq(T)
A52=3.6522e-5;
A53=0.16626;
A54=-4.56877*10;
A55=2.8733e4;
mu=A52+A53./T+A54./(T.^2)+A55./(T.^3);