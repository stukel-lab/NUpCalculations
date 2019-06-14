function [rhoreg,conf_rhoreg] = ...
    ComputeRhoRegMC(T,P,IpT,Ip0,Ispk,Iamb,Nspk,Namb,a,sig_T,sig_P,sig_IpT,sig_Ip0,sig_Ispk,sig_Iamb,sig_Nspk,sig_Namb,sig_a)

% This code calculates nutrient uptake rate (rhoreg), assuming nutrient 
% regeneration within the incubation bottle, using Eq. 13 in Stukel et 
% al. (submitted).  This code also calculates 
% asymmetric confidence limits for rhokan using a Monte Carlo approach.   

% Input parameters are:
% T = duration of incubation
% P = particulate concentration at end of incubation
% IpT = Isotope ratio of particulate pool at end of incubation
% Ip0 = Isotope ratio of particulate pool at beginning of incubation
% Ispk = Isotope ratio of tracer spike
% Iamb = Isotope ratio of ambient nutrient pool
% Nspk = Concentration of tracer spike
% Namb = Concentration of ambient nutrient pool
% a = ratio of nutrient regeneration to nutrient uptake
% sig_T = uncertainty in duration of incubation
% sig_P = uncertainty in particulate concentration at end of incubation
% sig_IpT = uncertainty in Isotope ratio of particulate pool at end of incubation
% sig_Ip0 = uncertainty in Isotope ratio of particulate pool at beginning of incubation
% sig_Ispk = uncertainty in Isotope ratio of tracer spike
% sig_Iamb = uncertainty in Isotope ratio of ambient nutrient pool
% sig_Nspk = uncertainty in Concentration of tracer spike
% sig_Namb = uncertainty in Concentration of ambient nutrient pool
% sig_a = uncertainty in ratio of nutrient regeneration to nutrient uptake


[rhoreg] = ComputeRhoReg(T,P,IpT,Ip0,Ispk,Iamb,Nspk,Namb,a,sig_T,sig_P,sig_IpT,sig_Ip0,sig_Ispk,sig_Iamb,sig_Nspk,sig_Namb,sig_a);

N = 1000; %Number of simulations to conduct

for i=1:N
    T2 = T + randn*sig_T;
    P2 = P + randn*sig_P;
    IpT2 = IpT + randn*sig_IpT;
    Ip02 = Ip0 + randn*sig_Ip0;
    Ispk2 = Ispk + randn*sig_Ispk;
    Iamb2 = Iamb + randn*sig_Iamb;
    Nspk2 = Nspk + randn*sig_Nspk;
    Namb2 = Namb + randn*sig_Namb;
    a2 = a + randn*sig_a;
    Is0 = (Ispk2.*Nspk2+Iamb2.*Namb2)./(Nspk2+Namb2);
    if T2<=0 | Ispk2<=Iamb2
        rho_temp(i)=Inf;
    elseif P2<=0 | IpT2<=Ip02 | Namb2<=0
        rho_temp(i)=0;
    elseif (Ip02-IpT2)/(Nspk2+Namb2)+(Is0-a2*IpT2)/P2 <=0
        rho_temp(i)=Inf;
    else
        temp=ComputeRhoReg(T2,P2,IpT2,Ip02,Ispk2,Iamb2,Nspk2,Namb2,a2,sig_T,sig_P,sig_IpT,sig_Ip0,sig_Ispk,sig_Iamb,sig_Nspk,sig_Namb,sig_a);
        rho_temp(i)=temp;
    end
end
rho_temp=sort(rho_temp);
conf_rhoreg = [rho_temp(round(N*0.15865)) rho_temp(round(N*0.84135))];