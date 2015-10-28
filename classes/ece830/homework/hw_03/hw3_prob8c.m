clear;
close all;

x = linspace(0, 20, 1000);

h0 = exp(-x);
h1 = 0.5 * exp(-x / 2.0);

%disp('Mean H0 =' + mean(h0));
%disp('Mean H1 =' + mean(h1));

gamma = -log(0.05);
gamma_1 = exp(gamma / 2.0);
gamma_2 = (1 - exp(-gamma / 2.0)) / (1 - exp(-gamma));

disp('gamma, gamma1, gamma2' )
disp(gamma)
disp(gamma_1)
disp(gamma_2)

P_D = exp(-gamma / 2.0)
P_FA = 0.05

E_0 = ((1 - P_FA) * log((1 - P_FA) / (1 - P_D)) - ...
       P_FA * log((P_D) / (P_FA))) / KLDiv(h0, h1)

E_1 = (P_D * log(P_D / P_FA) - ...
       (1 - P_D) * log((1 - P_FA) / (1 - P_D))) / KLDiv(h0, h1)




