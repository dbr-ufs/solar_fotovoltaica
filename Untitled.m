curvasIV = dlmread('curvasIV_T_cte.csv',',');
Vm = curvasIV(:,1);
Im = curvasIV(:,2);
Tm = curvasIV(:,3);
Gm = curvasIV(:,4);

figure
plot(Vm,Im, '.')
xlim([0 max(Vm)])
ylim([0 max(Im)])
title('dados de entrada');

T0 = 25; %temperatura stc
E0 = 1000;
E = Gm(find(~Im))';
Temp = Tm(find(~Im))';
%E = unique(Gm');
%Temp = unique(Tm');
k = 1.38066E-23; % J/K
q = 1.60218E-19; % c
Ns = 36; %número de células
Vth = [Ns*k*(Temp + 273.15)/q];
n = 1.2; % for cSi cell
Voc = Vm(find(~Im))';
Isc = Im(find(~Vm))';

%% Regressao para obter beta e Voc0
x = [ones(length(Temp),1), [Temp.' - T0]];
y = [Voc - n*Vth.*log(E/E0)].';

coeficiente = x\y;
Voc_regress = coeficiente(1);
beta = coeficiente(2);
beta_relativo = coeficiente(2)/Voc_regress*100;
fprintf('Resultados da regressao:\n');
fprintf('Voc_stc:       %.5f   V\n', Voc_regress);
fprintf('beta:          %.5f   V/ºC\n', beta);
fprintf('beta relativo: %.5f   %%/ºC\n\n',beta_relativo) 

xi = [ones(length(Temp),1), [Temp.' - T0]];
yi = [Isc.*(E0./E)];

coeficiente = xi\yi.';
Isc_regress = coeficiente(1);
alfa = coeficiente(2);
alfa_relativo = coeficiente(2)/Isc_regress*100;
fprintf('Resultados da regressao:\n');
fprintf('Isc_stc:       %.5f   A\n', Isc_regress);
fprintf('alfa:          %.5f   A/ºC\n', alfa);
fprintf('alfa relativo: %.5f   %%/ºC\n',alfa_relativo) 

% T e G dos valores máximos de irradiância
[Gobj,k] = max(Gm)
Tobj = Tm(k)

%% Determinacao de Kappa
Rs = 0; a = 0;
kappa_objfun = @(kappa) max_P_diff_for_Rs(Im,Vm,Tm,Gm, Tobj, Gobj, alfa, beta, Rs, a, kappa);
options = optimset('TolX', 1e-6);
kappa = fminsearch(kappa_objfun, 0, options);
fprintf('kappa = %.7f\n', kappa);
erro_kappa = kappa_objfun(kappa);
fprintf('erro_kappa  = %.7f\n', erro_kappa);
%kappa = 0.0032334; % [V/A.ºC]


%% correcao para Ponto de Máxima potência
function val = max_P_diff_for_Rs(Im,Vm,Tm,Gm, Tobj, Gobj, alfa, beta, Rs, a,kappa)
    Voc = Vm(find(~Im));
    Voc2 = Voc(1)*ones(56,1); % são 56 dados + Isc + Voc
    for i=2:size(Voc,1)
        Voc2 = [Voc2; Voc(i)*ones(56,1)];
    end
    Itransl = Im .* (1 + alfa*(Tobj - Tm)) .* Gobj ./ Gm;
    Vtransl = Vm + Voc2.*(beta.*(Tobj - Tm) + a.*log(Gobj./Gm)) - Rs.*(Itransl - Im) - kappa.*Itransl.*(Tobj - Tm);
    Pm = Vtransl .* Itransl;
    k = find(~Im);
    Pmp = max(Pm(1:k(1)));
    for i = 2:(size(k,1)-1)
        Pmp(i) = max(Pm(k(i):k(i+1)));
    end
    val = max(abs((Pmp - median(Pmp))./median(Pmp)));
end