clear
close all

m = 5;
n = 9;

a = rand(m, 1);
epsilons = [0.0001, 0.001, 0.01];

fig = figure(1);clf;
subplot(111); scatter(a, 1.0+2.0.*a); hold on;

b_fits = zeros(3,1);

for i=1:length(epsilons)
    epsilon = epsilons(i);

    I = eye(n);
    zero_matrix = zeros(n, 1);

    A = [[ones(m, 1), a, a.^2, a.^3, a.^4, a.^5, a.^6, a.^7, a.^8]; ...
         [epsilon.* I]];

    b = [[1.0+2.0.*a]; [zero_matrix]];

    x = A\b;

    a_fit = linspace(0,1,100);
    b_fit = x(1) + x(2).*a_fit + x(3).*a_fit.^2 + x(4).*a_fit.^3 + ...
            x(5).*a_fit.^4 + x(6).*a_fit.^5 ...
            + x(7).*a_fit.^7 + x(8).*a_fit.^8;

    plot(a_fit, b_fit); 
end
xlabel('a')
ylabel('b')
hold off;
%legend('data', 'eps = 0.001', 'eps = 0.01', 'eps = 0.1', 'Location', 'North')
box on;

saveas(fig, 'problem3c_fig', 'png')


