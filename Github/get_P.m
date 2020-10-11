function P=get_P(t,tc,Td,Tw,Tc)
S1=get_S(t,tc-Td/2,0,Tw,Tc);
S2=get_S(t,tc-Td/2,Td,Tw,Tc);
S3=get_S(t,tc-Td/2,Tc,Tw,Tc);
S4=get_S(t,tc-Td/2,Tc+Td,Tw,Tc);
P=S1-S2+S3-S4;
