close all;
clear all;

par.Gd=18e3;
par.T0=318;

par.Gv=60e3;
par.Cp=4d3;
par.TH=278;
par.Ts=283;
S=10;
temp=(260:350);
[fv,fE,fK]=calcfTs(temp,par);
kT=calcV(temp,par,S);

dlnk=calcdlnV(temp,par,S);
par.Gv=50e3;
kT1=calcV(temp,par,S);

par.Cp=10d3;
kT2=calcV(temp,par,S);

subplot(2,1,1);plot(temp,log(kT));
hold on;
plot(temp,log(kT1));
plot(temp,log(kT2));

subplot(2,1,2);plot(temp,(kT));
hold on;
plot(temp,kT1);
plot(temp,(kT2));

dlnk1=calcdlnV(temp,par,S);

