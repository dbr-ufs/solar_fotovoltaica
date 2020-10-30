curvasIV = dlmread('curvasIV.csv',',');
Vm = curvasIV(:,1);
Im = curvasIV(:,2);
Tm = curvasIV(:,3);
Gm = curvasIV(:,4);

figure
plot(Vm,Im, '*')
xlim([0 max(Vm)])
ylim([0 max(Im)])
title('dados de entrada');

% Coeficientes Correção
alfa =   0.06020000/100; % [1/ºC] na Corrente
beta =  -0.35468047/100; % [1/ºC] na Tensão
a  = 0.0647002;
Rs = 0.2779663; % [Ohm]
kappa = 0.0032334; % [V/A.ºC]

Gobj = 1000; % [W/m2]
Tobj = 25; % ºC
    
Itransl = Im .* (1 + alfa.*(Tobj - Tm)) .* Gobj ./ Gm;
Vtransl = Vm + max(Vm).*(beta.*(Tobj - Tm) + a.*log(Gobj./Gm)) - Rs.*(Itransl - Im) - kappa.*Itransl.*(Tobj - Tm);

figure
plot(Vtransl,Itransl, '.')
title('Todas as curvas corrigidas para stc')
