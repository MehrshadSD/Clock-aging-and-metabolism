function S=get_S(t,tf,ts,Tw,Tc)
tau=get_tau(t-tf,Tc);
S=1/2*(1+tanh((tau-ts)/Tw));
