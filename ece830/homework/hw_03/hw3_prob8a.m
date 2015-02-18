clear;
close all;

x = linspace(0, 20, 1000);

h0 = exp(-x);
h1 = 0.5 * exp(-x / 2.0);

mean(h0)
mean(h1)

% Plot PDFs
figure(1); 
subplot(121); hold on;
plot(x, h0, 'k.');
axis('square');
title('H0');
xlabel('x');
ylabel('PDF');

subplot(122); hold on;
plot(x, h1, 'k.');
axis('square');
title('H0');
xlabel('x');

hold off;

