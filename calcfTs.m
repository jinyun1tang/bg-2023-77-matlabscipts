function [fv,fE,fK]=calcfTs(temp,par)
%compute the temperature factors

RTi=1./(8.314.*temp);
fK=exp(-par.GK.*RTi.*(1.-temp./par.T0));
fv=temp./par.T0.*exp(-par.Gv.*RTi.*(1.-temp./par.T0));
DGE=par.Cp.*((temp-par.TH)+temp.*log(par.Ts./temp));
fE=1./(1.+exp(-DGE.*RTi));


end