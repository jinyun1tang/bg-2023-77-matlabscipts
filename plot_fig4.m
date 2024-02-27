close all;
clear all;
disp_flag=0;
example=1;
temp=(250:380);   
Extra.S=0.5;

[par,obs,stat,model]=fParEstimate(example,disp_flag);
[fv,fE,fK]=calcfTs(temp,par);
vmax=fv.*fE;
T3=temp(vmax==max(vmax));    

par.Gv=par.Gv-20e3;
[fv,fE,fK]=calcfTs(temp,par);
vmax1=fv.*fE;
maxv=max(vmax).*0.01;
T4=temp(vmax1==max(vmax1));
plot(temp,vmax./maxv,'LineWidth',2);
hold on;
plot(temp,vmax1./maxv,'LineWidth',2);
set(gcf,'color','w');
set(gca,'FontSize',14);
xlabel('Temperature (K)','FontSize',18);
ylabel('% of maximum reaction rate','FontSize',18);
ch=legend('High activation energy: T_o_p_t=330 K','Low activation energy: T_o_p_t=328 K');
set(ch,'FontSize',16,'location','northwest');