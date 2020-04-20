addpath functions

load('data/RCWA_data.mat', 'kx', 'lambda', 'data')

h = 1.62; % slab thickness
bw_cutoff = 3; % see BICModel help
max_p = 2; % max phase-matching order
[ld0, Q] = BICModel(kx, lambda, data, h, bw_cutoff, max_p);

I = find(Q{2}>3e7);
figure
subplot(1,2,1)
hold on
plot(kx/2/pi, ld0{2}, 'g-', 'LineWidth', 2)
plot(kx(I)/2/pi, ld0{2}(I), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g')
ylabel 'a/\lambda'
xlabel 'k_xa/(2\pi)'

subplot(1,2,2)
semilogy(kx/2/pi, Q{2}, 'g-', 'LineWidth', 2)
ylabel Q
xlabel 'k_xa/(2\pi)'

set(get(gcf, 'Children'), 'FontSize', 14);
