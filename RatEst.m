close all;
clear all;
clc;
%Cp, H*
Rat=[75.54 5341;
81.36 5311;
78.48 5331;
75.1 5354;
77.27 5338;
80.92 5320;
75.78 5342;
75.38 5344;
73.31 5347;
79.68 5329;
74.56 5353;
73.9 5353;
74.59 5353;
73.19 5351;
73.59 5354;
73.68 5355;
75.2 5351;
75.43 5348;
73.64 5354;
71.87 5360;
70.61 5362;
71.09 5354;
70.87 5354;
69.1 5360;
68.98 5359;
70.79 5354;
71.09 5362;
64.54 5365;
63 5355;
63.06 5358;
57.76 5324;
58 5327;
57.54 5326;
63.41 5338;
56.58 5318];

THs=373.6;
TH=-Rat(:,2)./Rat(:,1)+THs;
Tss=385.2;
TS=Tss.*exp(-18.1./Rat(:,1));

for jj=1:size(Rat,1)
    fprintf('%d: %f,%f,%f\n',jj,TH(jj),TS(jj),TS(jj)-TH(jj));
end
    