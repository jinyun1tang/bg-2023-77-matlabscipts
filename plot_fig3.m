close all;
clear all;
addpath('/Users/jinyuntang/work/github/matlab_tools');
disp_flag=0;
fig=figure(1);
ax=multipanel(fig,3,3,[0.06,0.08],[0.26,0.27],[0.05,0.05]);
temp=(250:380);   
Extra.S=0.5;
example1=0;
labels={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) '};
for example=1 : 12
    if((example-1)*(example-7)*(example-9)*(example-11)==0)
        continue
    end
    example1=example1+1
    [par,obs,stat,model]=fParEstimate(example,disp_flag);

    kT0=calcV(temp,par,Extra.S);
    kT1=calcV(temp,par,Extra.S.*2);
    kT2=calcV(temp,par,Extra.S.*20);    
    [fv,fE,fK]=calcfTs(temp,par);
    vmax=fv.*fE;

    T0=temp(kT0==max(kT0));
    T1=temp(kT1==max(kT1));
    T2=temp(kT2==max(kT2));
    T3=temp(vmax==max(vmax));    

    set_curAX(fig,ax(example1));
    plot(temp,kT0./max(vmax).*100,'LineWidth',2);
    hold on;
    plot(temp,kT1./max(vmax).*100,'LineWidth',2);
    plot(temp,kT2./max(vmax).*100,'LineWidth',2);
    plot(temp,vmax./max(vmax).*100,'LineWidth',2);

    T0s=sprintf('T_o_p_t_,_c_1=%.0f K',T0);
    T1s=sprintf('T_o_p_t_,_c_2=%.0f K',T1);
    T2s=sprintf('T_o_p_t_,_c_3=%.0f K',T2);    
    T3s=sprintf('T_o_p_t_,_c_4=%.0f K',T3);
    ylim([0,120]);
    put_tag(fig,ax(example1),[0.1,.90],[labels{example1}, stat.name],16);
    put_tag(fig,ax(example1),[0.1,.78],T0s,16);
    put_tag(fig,ax(example1),[0.1,.66],T1s,16);
    put_tag(fig,ax(example1),[0.1,.54],T2s,16);
    put_tag(fig,ax(example1),[0.1,.42],T3s,16);
    set(ax,'FontSize',16,'XTick',(250:50:400));
    if(example1>=6)
       xlabel('Temperature (K)','FontSize',16);        
    end
    if(mod(example1,3)==1)
       ylabel('% of maximum reaction rate','FontSize',16);
    end

end
ch=legend('c1) S=K_0/2','c2) S=K_0','c3) S=10K_0','c4) S=\infty');
delete(ax(9));

fig=figure(2);
ax=multipanel(fig,2,2,[0.06,0.08],[0.4,0.4],[0.05,0.075]);

example=1;
[par,obs,stat,model]=fParEstimate(example,disp_flag);
kT0=calcV(temp,par,Extra.S);
kT1=calcV(temp,par,Extra.S.*2);
kT2=calcV(temp,par,Extra.S.*20);    
[fv,fE,fK]=calcfTs(temp,par);
vmax=fv.*fE;
T0=temp(kT0==max(kT0));
T1=temp(kT1==max(kT1));
T2=temp(kT2==max(kT2));
T3=temp(vmax==max(vmax));    

set_curAX(fig,ax(1));
h(1)=plot(temp,kT0./max(vmax).*100,'LineWidth',2);
hold on;
h(2)=plot(temp,kT1./max(vmax).*100,'LineWidth',2);
h(3)=plot(temp,kT2./max(vmax).*100,'LineWidth',2);
h(4)=plot(temp,vmax./max(vmax).*100,'LineWidth',2);

ch=legend(h,'c1) S=K_0/2','c2) S=K_0','c3) S=10K_0','c4) S=\infty');
set(ch,'FontSize',16);
ylim([0,120]);

T0s=sprintf('T_o_p_t_,_c_1=%.0f K',T0);
T1s=sprintf('T_o_p_t_,_c_2=%.0f K',T1);
T2s=sprintf('T_o_p_t_,_c_3=%.0f K',T2);
T3s=sprintf('T_o_p_t_,_c_4=%.0f K',T3);

put_tag(fig,ax(1),[0.05,.9],['(a) ', stat.name],16);
put_tag(fig,ax(1),[0.05,.78],T0s,16);
put_tag(fig,ax(1),[0.05,.66],T1s,16);
put_tag(fig,ax(1),[0.05,.54],T2s,16);
put_tag(fig,ax(1),[0.05,.42],T3s,16);

example=7;
[par,obs,stat,model]=fParEstimate(example,disp_flag);
kT0=calcV(temp,par,Extra.S);
kT1=calcV(temp,par,Extra.S.*2);
kT2=calcV(temp,par,Extra.S.*20);    
[fv,fE,fK]=calcfTs(temp,par);
vmax=fv.*fE;
T0=temp(kT0==max(kT0));
T1=temp(kT1==max(kT1));
T2=temp(kT2==max(kT2));
T3=temp(vmax==max(vmax));       
set_curAX(fig,ax(2));
plot(temp,kT0./max(vmax).*100,'LineWidth',2);
hold on;
plot(temp,kT1./max(vmax).*100,'LineWidth',2);
plot(temp,kT2./max(vmax).*100,'LineWidth',2);
plot(temp,vmax./max(vmax).*100,'LineWidth',2);

ylim([0,120]);

T0s=sprintf('T_o_p_t_,_c_1=%.0f K',T0);
T1s=sprintf('T_o_p_t_,_c_2=%.0f K',T1);
T2s=sprintf('T_o_p_t_,_c_3=%.0f K',T2);
T3s=sprintf('T_o_p_t_,_c_4=%.0f K',T3);

put_tag(fig,ax(2),[0.05,.9],['(b) ', stat.name],16);
put_tag(fig,ax(2),[0.05,.78],T0s,16);
put_tag(fig,ax(2),[0.05,.66],T1s,16);
put_tag(fig,ax(2),[0.05,.54],T2s,16);
put_tag(fig,ax(2),[0.05,.42],T3s,16);

example=9;
[par,obs,stat,model]=fParEstimate(example,disp_flag);
kT0=calcV(temp,par,Extra.S);
kT1=calcV(temp,par,Extra.S.*2);
kT2=calcV(temp,par,Extra.S.*20);    
[fv,fE,fK]=calcfTs(temp,par);
vmax=fv.*fE;
T0=temp(kT0==max(kT0));
T1=temp(kT1==max(kT1));
T2=temp(kT2==max(kT2));
T3=temp(vmax==max(vmax));       
    
set_curAX(fig,ax(3));
plot(temp,kT0./max(vmax).*100,'LineWidth',2);
hold on;
plot(temp,kT1./max(vmax).*100,'LineWidth',2);
plot(temp,kT2./max(vmax).*100,'LineWidth',2);
plot(temp,vmax./max(vmax).*100,'LineWidth',2);

ylim([0,120]);

T0s=sprintf('T_o_p_t_,_c_1=%.0f K',T0);
T1s=sprintf('T_o_p_t_,_c_2=%.0f K',T1);
T2s=sprintf('T_o_p_t_,_c_3=%.0f K',T2);
T3s=sprintf('T_o_p_t_,_c_4=%.0f K',T3);

put_tag(fig,ax(3),[0.05,.9],['(c) ', stat.name],16);
put_tag(fig,ax(3),[0.05,.78],T0s,16);
put_tag(fig,ax(3),[0.05,.66],T1s,16);
put_tag(fig,ax(3),[0.05,.54],T2s,16);
put_tag(fig,ax(3),[0.05,.42],T3s,16);

example=11;
[par,obs,stat,model]=fParEstimate(example,disp_flag);
kT0=calcV(temp,par,Extra.S);
kT1=calcV(temp,par,Extra.S.*2);
kT2=calcV(temp,par,Extra.S.*20);    
[fv,fE,fK]=calcfTs(temp,par);
vmax=fv.*fE;
T0=temp(kT0==max(kT0));
T1=temp(kT1==max(kT1));
T2=temp(kT2==max(kT2));
T3=temp(vmax==max(vmax));       
    
set_curAX(fig,ax(4));
plot(temp,kT0./max(vmax).*100,'LineWidth',2);
hold on;
plot(temp,kT1./max(vmax).*100,'LineWidth',2);
plot(temp,kT2./max(vmax).*100,'LineWidth',2);
plot(temp,vmax./max(vmax).*100,'LineWidth',2);

ylim([0,120]);

T0s=sprintf('T_o_p_t_,_c_1=%.0f K',T0);
T1s=sprintf('T_o_p_t_,_c_2=%.0f K',T1);
T2s=sprintf('T_o_p_t_,_c_3=%.0f K',T2);
T3s=sprintf('T_o_p_t_,_c_4=%.0f K',T3);

put_tag(fig,ax(4),[0.05,.9],['(d) ', stat.name],16);
put_tag(fig,ax(4),[0.05,.78],T0s,16);
put_tag(fig,ax(4),[0.05,.66],T1s,16);
put_tag(fig,ax(4),[0.05,.54],T2s,16);
put_tag(fig,ax(4),[0.05,.42],T3s,16);

set(ax,'FontSize',16);

set_curAX(fig,ax(1));
ylabel('% of maximum reaction rate','FontSize',16);
set_curAX(fig,ax(3));
ylabel('% of maximum reaction rate','FontSize',16);
xlabel('Temperature (K)','FontSize',16);
set_curAX(fig,ax(4));
xlabel('Temperature (K)','FontSize',16);
