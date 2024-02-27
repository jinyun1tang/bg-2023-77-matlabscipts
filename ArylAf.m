function jcost=ArylAf(x,Extra)
%computing likelhiood log for the q-y function



par.Gv=x(1);
par.GK=x(1)-18d3;
par.Cp=x(2);
par.TH=x(3);
par.Ts=x(3)+x(4);
par.T0=Extra.T0;

dlnk=calcdlnV(Extra.ArylA(:,1),par,Extra.S);

jcost=sum((dlnk-log(Extra.ArylA(:,2))+log(Extra.ArylA(9,2))).^2,'omitnan')+exp(-10.*(par.Ts-par.TH));%+exp(-1.e-1.*(x(1)-x(2)))+exp(-0.1.*min(x(1:2)));
end