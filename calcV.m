function kT=calcV(temp,par,S)
%normalized velocity
[fv,fE,fK]=calcfTs(temp,par);
kT=fv.*fE.*S./(fK+S);
end