%Main
clear all; close all;
NPt = 5000;%10000;
for Adr=[0.028e-6];%0.028e-6]%0.1e-6,0.21e-6,0.5e-6,1e-6,2e-6,4e-6,5e-6];%0.5e-6,1e-6,2e-6,4e-6,5e-6
       Gross_Pitaevskii__Process(1e-7,Adr,NPt);
end
