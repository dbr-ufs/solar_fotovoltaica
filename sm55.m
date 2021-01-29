function Ia = sm55(Va,Rad,TaC)
% ms50.m model for the MS55 solar array
% current given voltage, illumination and temperature
% Ia = ms55(Va,G,TaC) = array voltage
% Ia,Va = array current,voltage
% G = num of Suns (1 Sun = 1000 W/m^2)
% TaC = Temp ambient in Deg C
if Rad > 5
    Suns = Rad / 1000;
else
    Suns = 0.005;
end
k = 1.38e-23; % Boltzman´s const
q = 1.60e-19; % charge on an electron
% enter the following constants here, and the model will be
% calculated based on these. for 1000W/m^2
A = 1.3;%1.04;
% "diode quality" factor, =2 for crystaline, <2 for amorphous
Vg = 1.12; % band gap voltage, 1.12eV for xtal Si, ~1.75 for amorphous Si.
Ns = 36; % number of series connected cells (diodes)
T1 = 273 + 25;
Voc_T1 = 21.7 /Ns; % open cct voltage per cell at temperature T1
Isc_T1 = 3.5; % short cct current per cell at temp T1
T2 = 273 + 60;
Voc_T2 = 18.30 /Ns; % open cct voltage per cell at temperature T2
Isc_T2 = 3.45; % short cct current per cell at temp T2
TaK = 273 + TaC;% + 35 * Suns; % array working temp
TrK = 273 + 25; % reference temp
% when Va = 0, light generated current Iph_T1 = array short cct current
% constant "a" can be determined from Isc vs T
Iph_T1 = Isc_T1 * Suns;
%a = (Isc_T2 - Isc_T1)/Isc_T1 * 1/(T2 - T1);
a = 0.0004; % Temperature Coefficients of Isc
Iph = Iph_T1 * (1 + a*(TaK - T1));
Vt_T1 = k * T1 / q; % = A * kT/q
Ir_T1 = Isc_T1 / (exp(Voc_T1/(A*Vt_T1))-1);
Ir_T2 = Isc_T2 / (exp(Voc_T2/(A*Vt_T1))-1);
b = Vg * q/(A*k);
Ir = Ir_T1 * (TaK/T1).^(3/A) .* exp(-b.*(1./TaK - 1/T1));
%Ir = 4.76e-10; %%%%%%%%%%%%%%%%

X2v = Ir_T1/(A*Vt_T1) * exp(Voc_T1/(A*Vt_T1));
dVdI_Voc = - 1.15/Ns / 2; % dV/dI at Voc per cell --
% from manufacturers graph
Rs = - dVdI_Voc - 1/X2v; % series resistance per cell
% Ia = 0:0.01:Iph;
Vt_Ta = A * k * TaK / q; % = A * kT/q
% Ia1 = Iph - Ir.*( exp((Vc+Ia.*Rs)./Vt_Ta) -1);
% solve for Ia: f(Ia) = Iph - Ia - Ir.*( exp((Vc+Ia.*Rs)./Vt_Ta) -1) = 0;
% Newton´s method: Ia2 = Ia1 - f(Ia1)/f´(Ia1)
Vc = Va/Ns;
Ia = zeros(size(Vc));
Rs = 0.476/ Ns;
Rp = 141/ Ns;
%Iav = Ia;
for j=1:5;
Ia = Ia - ...
(Iph - Ia - Ir.*( exp((Vc+Ia.*Rs)./Vt_Ta) -1) - (Vc+Ia.*Rs)./Rp...
)./ (-1 - (Ir.*( exp((Vc+Ia.*Rs)./Vt_Ta) -1)).*Rs./Vt_Ta - Rs./Rp);
 %Iav = [Iav;Ia]; % to observe convergence for debugging.
end
betaRef = -0.0045;
%Rs = 0.476/ Ns;
T_cell = TaK;
Voc_STC = Voc_T1;
Tref = T1;
%A = 1;
%Vc = Vc * Ns;
for j=1:50
    Vt_Ta = A * k * T_cell / q;
    fT = (Ia + Ir.*( exp((Vc+Ia.*Rs)./Vt_Ta) -1) + (Vc+Ia.*Rs)./Rp - ...
        Ir.*( exp((Voc_STC*(1 + betaRef*(T_cell - Tref)) + A*k/q*T_cell*log(Suns))./Vt_Ta) -1) ...
        ).* Rp - (Voc_STC*(1 + betaRef*(T_cell - Tref)) + A*k/q*T_cell*log(Suns));
    f_T = (Ir.* exp(Voc_STC*(1 + betaRef*(T_cell - Tref))./Vt_Ta) ...
        .* ( Voc_STC*(1 + betaRef*(T_cell - Tref))./Vt_Ta./T_cell ...
        - (Voc_STC *betaRef - A*k/q*log(Suns)) ./ Vt_Ta) - ...
        Ir.* exp((Vc+Ia.*Rs)./Vt_Ta).* (Vc+Ia.*Rs)./Vt_Ta./T_cell ...
        ).* Rp - Voc_STC * betaRef - A*k/q*log(Suns);
    T_cell = T_cell - fT ./ f_T;
end
T_cell-273.15