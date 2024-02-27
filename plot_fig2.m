close all;
clear all;
addpath('/Users/jinyuntang/work/github/matlab_tools');
format long;
fig=figure(1);
ax=multipanel(fig,3,4,[0.06,0.08],[0.18,0.265],[0.05,0.05]);
tags={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'};
disp_flag=0;
for example=1 : 12
[par,obs,stat,model]=fParEstimate(example,disp_flag);
set_curAX(fig,ax(example));
plot(obs(:,1),obs(:,2),'o','MarkerFaceColor','r','MarkerSize',8);
hold on;
plot(model(:,1),model(:,2),'b-','LineWidth',2);

if mod(example,4)==1
    ylabel('ln F(T)- ln F(T_0)','FontSize',16);
end
if example > 8
    xlabel('Temperature (K)','FontSize',16);
end
ylim([-6,3]);
put_tag(fig,ax(example),[0.015,0.95],[tags{example},' ',stat.name],16);
Gv=sprintf('\\DeltaH_v=%.2f kJ mol^-^1',par.Gv.*1.e-3);
Cp=sprintf('\\DeltaC_p=%.2f kJ mol^-^1 K^-^1',par.Cp*1.e-3);
Ts=sprintf('T_S=%.1f K',par.Ts);
TH=sprintf('T_H=%.1f K',par.TH);
T0=sprintf('T_r=%.1f K',par.T0);
R2=sprintf('R^2=%.2f',stat.R2);
put_tag(fig,ax(example),[0.38,0.46],TH,16);
put_tag(fig,ax(example),[0.38,0.33],Ts,16);
put_tag(fig,ax(example),[0.17,0.20],Gv,16);
put_tag(fig,ax(example),[0.17,0.09],Cp,16);
put_tag(fig,ax(example),[0.09,0.81],T0,16);
put_tag(fig,ax(example),[0.75,0.92],R2,16);

end

set(ax,'FontSize',16);

