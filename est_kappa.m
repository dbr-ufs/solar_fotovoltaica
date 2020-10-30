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
%kappa = 0.0032334; % [V/A.ºC]

Gobj = 1000; % [W/m2]
Tobj = 25; % ºC

%% Determinacao de Kappa
%Rs = 0; a = 0;
kappa_objfun = @(kappa) max_P_diff_for_Rs(Im,Vm,Tm,Gm, Tobj, Gobj, alfa, beta, Rs, a, kappa);
options = optimset('TolX', 1e-6);
kappa = fminsearch(kappa_objfun, 0, options);
fprintf('kappa = %.7f\n', kappa);
erro_kappa = kappa_objfun(kappa);
fprintf('erro_kappa  = %.7f\n', erro_kappa);
%% Resultado

Itransl = Im .* (1 + alfa.*(Tobj - Tm)) .* Gobj ./ Gm;
Vtransl = Vm + max(Vm).*(beta.*(Tobj - Tm) + a.*log(Gobj./Gm)) - Rs.*(Itransl - Im) - kappa.*Itransl.*(Tobj - Tm);

figure
plot(Vtransl,Itransl, '.')
title('Todas as curvas corrigidas para stc')


%% correcao para Ponto de Máxima potência
function val = max_P_diff_for_Rs(Im,Vm,Tm,Gm, Tobj, Gobj, alfa, beta, Rs, a,kappa)
    Itransl = Im .* (1 + alfa.*(Tobj - Tm)) .* Gobj ./ Gm;
    Vtransl = Vm + max(Vm).*(beta.*(Tobj - Tm) + a.*log(Gobj./Gm)) - Rs.*(Itransl - Im) - kappa.*Itransl.*(Tobj - Tm);
    for i=1:3 %4
        Pmp(i)     = max(Vtransl((1+56*(i-1)):56*i) .* Itransl((1+56*(i-1)):56*i));
    end
    val = max(abs((Pmp - median(Pmp))./median(Pmp)));
end

