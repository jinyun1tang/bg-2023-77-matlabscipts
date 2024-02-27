function dlnk=calcdlnV(temp,par,S)
%normalized velocity
[fv,fE,fK]=calcfTs(temp,par);
[fv0,fE0,fK0]=calcfTs(par.T0,par);

dlnk=log(fv.*fE./(fK+S))-log(fv0.*fE0./(fK0+S));
end