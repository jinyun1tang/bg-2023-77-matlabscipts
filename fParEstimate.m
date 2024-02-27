function [par,obs,stat,model]=fParEstimate(example,disp_flag)

Extra.opt=0;
if example==1
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/KT_data/V200S.mat');
    Extra.V200S=V200S;
    Extra.S=10;
    Extra.T0=V200S(9,1);
    fun = @(x)V200Sf(x,Extra);
elseif example==2
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/KT_data/V200A.mat');
    Extra.V200A=V200A;
    Extra.S=10;
    Extra.T0=V200A(9,1);
    fun = @(x)V200Af(x,Extra);    
elseif example==3
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/KT_data/G202P.mat');
    Extra.G202P=G202P;
    Extra.S=10;
    Extra.T0=G202P(3,1);
    fun = @(x)kG202Pf(x,Extra);    
elseif example==4
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/KT_data/WMaL.mat');
    Extra.WMaL=WMaL;
    Extra.S=10;
    Extra.T0=WMaL(8,1);
    fun = @(x)WMaLf(x,Extra);    
elseif example==5
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/KT_data/V200T.mat');
    Extra.V200T=V200T;
    Extra.S=10;
    Extra.T0=V200T(9,1);
    fun = @(x)V200Tf(x,Extra);    
elseif example==6
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/KT_data/WldB.mat');
    Extra.WldB=WldB;
    Extra.S=2;
    Extra.T0=WldB(5,1);
    fun = @(x)WldBf(x,Extra);       
elseif example==7
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/KT_data/A43b.mat');
    Extra.A43b=A43b;
    Extra.S=2;
    Extra.T0=A43b(7,1);
    fun = @(x)A43bf(x,Extra);            
elseif example==8
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/peterson_data/Acdp.mat');
    Extra.Acdp=Acdp;
    Extra.S=10;
    Extra.T0=Acdp(4,1);
    fun = @(x)Acdpf(x,Extra);            
elseif example==9
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/peterson_data/Aded.mat');
    Extra.Aded=Aded;
    Extra.S=10;
    Extra.T0=Aded(3,1);
    fun = @(x)Adedf(x,Extra);            
elseif example==10
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/peterson_data/AlkP.mat');
    Extra.Alkp=AlkP;
    Extra.S=10;
    Extra.T0=AlkP(5,1);
    fun = @(x)Alkpf(x,Extra);            
elseif example==11
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/peterson_data/ArylA.mat');
    Extra.ArylA=ArylA;
    Extra.S=10;
    Extra.T0=ArylA(9,1);
    fun = @(x)ArylAf(x,Extra);            
elseif example==12
    load('/Users/jinyuntang/work/github/paper_matlabscripts/2023/MMRT_foundation/peterson_data/BetL.mat');
    Extra.BetL=BetL;
    Extra.S=10;
    Extra.T0=BetL(3,1);
    fun = @(x)BetLf(x,Extra);            
end
x0=[60d3,10d3,272,5];

options = optimset('TolX',1e-7);
x = fminsearch(fun,x0,options);

hmat=hessianmat(fun,x)
inv(hmat)
switch example
    case 1
        stat.name='V200S';
        obs=[V200S(:,1),V200S(:,2)-V200S(9,2)];
    case 2
        stat.name='V200A';
        obs=[V200A(:,1),V200A(:,2)-V200A(9,2)];
    case 3    
        stat.name='G202P';
        obs=[G202P(:,1),G202P(:,2)-G202P(3,2)];         
    case 4    
        stat.name='Wildtype MalL';
        obs=[WMaL(:,1),WMaL(:,2)-WMaL(8,2)];    
    case 5
        stat.name='V200T';
        obs=[V200T(:,1),V200T(:,2)-V200T(9,2)];        
    case 6    
        stat.name='Barnase';
        obs=[WldB(:,1),WldB(:,2)-WldB(5,2)];    
    case 7    
        stat.name='A43C/S80C';
        obs=[A43b(:,1),A43b(:,2)-A43b(7,2)];     
    case 8    
        stat.name='Acid phosphatase';
        obs=[Acdp(:,1),log(Acdp(:,2))-log(Acdp(4,2))];         
    case 9    
        stat.name='Adenosine deaminase';
        obs=[Aded(:,1),log(Aded(:,2))-log(Aded(3,2))];         
    case 10    
        stat.name='Alkaline phosphatase';
        obs=[AlkP(:,1),log(AlkP(:,2))-log(AlkP(5,2))];         
    case 11    
        stat.name='Aryl-acylamidase';
        obs=[ArylA(:,1),log(ArylA(:,2))-log(ArylA(9,2))];         
    case 12    
        stat.name='\beta-lactamse';
        obs=[BetL(:,1),log(BetL(:,2))-log(BetL(3,2))];         
end
par.GK=x(1)-18d3;
par.Gv=x(1);
par.Cp=x(2);
par.TH=x(3);
par.Ts=x(3)+x(4);
par.T0=Extra.T0;
temp=(280:360);

dlnk=calcdlnV(temp,par,Extra.S);
model=[temp',dlnk'];
dlnk1=calcdlnV(obs(:,1),par,Extra.S);
w = linearfit(dlnk1,obs(:,2),0.05,'disp');
stat.R2=w(6);
stat.RMSE=w(5);
if(disp_flag)
    subplot(2,2,1);
    plot(temp,dlnk);
    hold on;
   
   
if example==1    
    plot(V200S(:,1),V200S(:,2)-V200S(9,2),'o','MarkerFaceColor','r');
elseif example==2
    plot(V200A(:,1),V200A(:,2)-V200A(9,2),'o','MarkerFaceColor','r');    
elseif example==3
    plot(G202P(:,1),G202P(:,2)-G202P(3,2),'o','MarkerFaceColor','r');    
elseif example==4
    plot(WMaL(:,1),WMaL(:,2)-WMaL(8,2),'o','MarkerFaceColor','r');    
elseif example==5
    plot(V200T(:,1),V200T(:,2)-V200T(9,2),'o','MarkerFaceColor','r');        
elseif example==6
    plot(WldB(:,1),WldB(:,2)-WldB(5,2),'o','MarkerFaceColor','r');    
elseif example==7
    plot(A43b(:,1),A43b(:,2)-A43b(7,2),'o','MarkerFaceColor','r');     
elseif example==8
    plot(Acdp(:,1),log(Acdp(:,2))-log(Acdp(4,2)),'o','MarkerFaceColor','r');         
elseif example==9
    plot(Aded(:,1),log(Aded(:,2))-log(Aded(3,2)),'o','MarkerFaceColor','r');         
elseif example==10
    plot(AlkP(:,1),log(AlkP(:,2))-log(AlkP(5,2)),'o','MarkerFaceColor','r');         
elseif example==11
    plot(ArylA(:,1),log(ArylA(:,2))-log(ArylA(9,2)),'o','MarkerFaceColor','r');         
elseif example==12
    plot(BetL(:,1),log(BetL(:,2))-log(BetL(3,2)),'o','MarkerFaceColor','r');         
end
end
temp=(250:380);
kT=calcV(temp,par,Extra.S);
kT1=calcV(temp,par,Extra.S.*0.5);
[fv,fE,fK]=calcfTs(temp,par);

if(disp_flag)
subplot(2,2,2);plot(temp,kT./max(kT).*100);hold on;plot(temp,kT1./max(kT).*100);
if(example==8)
    plot(Acdp(:,1),Acdp(:,2),'o','MarkerFaceColor','r');
elseif(example==9)
    plot(Aded(:,1),Aded(:,2),'o','MarkerFaceColor','r');
elseif(example==10)
    plot(AlkP(:,1),AlkP(:,2),'o','MarkerFaceColor','r');
elseif(example==11)
    plot(ArylA(:,1),ArylA(:,2),'o','MarkerFaceColor','r');
elseif(example==12)
    plot(BetL(:,1),BetL(:,2),'o','MarkerFaceColor','r');
end
subplot(2,2,3);plot(temp,fE);
subplot(2,2,4);plot(temp,fK);
end
end