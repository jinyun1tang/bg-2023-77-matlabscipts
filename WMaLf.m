function jcost=WMaLf(x,Extra)
%computing likelhiood log for the q-y function



par.Gv=x(1);
par.GK=x(1)-18d3;
par.Cp=x(2);
par.TH=x(3);
par.Ts=x(3)+x(4);
par.T0=Extra.T0;

dlnk=calcdlnV(Extra.WMaL(:,1),par,Extra.S);

jcost=sum((dlnk-Extra.WMaL(:,2)+Extra.WMaL(8,2)).^2,'omitnan')+exp(-10.*(par.Ts-par.TH));%+exp(-0.05.*(x(1)-x(2)))+exp(-0.1.*min(x(1:2)));
if Extra.opt==1
jcost=-jcost;
end
end