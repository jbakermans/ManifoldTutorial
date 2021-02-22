% PART 2 - https://www.crowdcast.io/e/sara-sollas-world-wide, 45:00-55:00
% Load synthetic data from two different groups of cells, 
% to simulate different recording days, then find how latent dynamics
% across days are oriented relative to each other.
% Jacob Bakermans, February 2021
close all;
clear all;
%% Section 1: Prepare data
% In this section we'll generate some synthetic neural data for two different
% recording days, with three different neurons recorded on each day

% Generate synthetic data for two different days 
[days{1}, days{2}] = generate_data_2(); % Data matrix: N_neurons x N_timebins
% Remove mean in both
for currDay = 1:length(days)
    days{currDay} = days{currDay} - ...
        repmat(mean(days{currDay}, 2), [1, size(days{currDay},2)]);
end

% Cells: assume some overlap, two cells shared across days
cells{1} = [1, 2, 3];
cells{2} = [1, 2, 4];

% Prepare data matrices for union of all neurons across days for each day
D = length(unique(cat(1, cells{:})));
T = size(days{1}, 2);
X_n = zeros(D, T);
X_m = zeros(D, T);
for currCell = 1:length(cells{1})
    X_n(cells{1}(currCell), :) = days{1}(currCell, :);
    X_m(cells{2}(currCell), :) = days{2}(currCell, :);
end

% Plot data of both days
figure();
for currDay = 1:2
    % First subplot: neuron firing rate through time
    subplot(2, 2, 1 + (currDay - 1) *2);
    hold on;
    for currNeuron = 1:length(cells{currDay})
        plot(days{currDay}(currNeuron, :));
    end
    hold off;
    % Set plot layout properties
    legend(['Neuron ' num2str(cells{currDay}(1))],...
        ['Neuron ' num2str(cells{currDay}(2))], ...
        ['Neuron ' num2str(cells{currDay}(3))]);
    ylabel('Normalised spike counts (1)');
    xlabel('Time (bins)');
    title(['Day ' num2str(currDay) ': Neuron activity through time']);
    % Second subplot: trajectory in neural space
    subplot(2, 2, 2 + (currDay - 1) *2);
    % Plot data
    plot3(days{currDay}(1, :), ...
        days{currDay}(2, :),...
        days{currDay}(3, :), 'x-');
    % Set plot layout properties
    xlabel(['Neuron ' num2str(cells{currDay}(1)) ' activity']);
    ylabel(['Neuron ' num2str(cells{currDay}(2)) ' activity']);
    zlabel(['Neuron ' num2str(cells{currDay}(3)) ' activity']);
    view(75, 30);
    grid on;
    title(['Trajectory in neural space of day ' num2str(currDay)]);
end

%% Section 2: singular value decomposions
% In this section, we'll use singular value decomposition to find the low-
% dimensional manifolds in the ambient dimension D for the union of neurons
% across days. We'll plot the manifolds in the original neural activity 
% spaces because they are 3D, while the full ambient space is 4D

% Do singular value decomposition for both days. U contains basis vectors in
% D dimension in its columns and provides a basis in the space of neurons.
[U_n, S_n, V_n] = svd(X_n);
[U_m, S_m, V_m] = svd(X_m);

% Now keep only a few of the basis vectors in U: the first d dimensions, 
% corresponding to leading eigenvalues (Matlab's svd function orders basis 
% vectors by eigenvalues automatically). d is the dimensionality 
% of the low-dimensional neural manifold embedded in emperical neural space.

% Set 'flat dimension' d
d = 2;

% Select first d basis vectors from U_n and U_m
U_n_tilde = U_n(:, 1:d);
U_m_tilde = U_m(:, 1:d);

% Plot matrices we've built so far
figure();
% First subplot: U_n
subplot(2,2,1);
imagesc(U_n);
axis image;
xlabel('Columns (basis vectors)');
ylabel('Rows (neural space dimensions)');
title('$U_n$', 'interpreter','latex');
% Second subplot: U_m
subplot(2,2,2);
imagesc(U_m);
axis image;
xlabel('Columns (basis vectors)');
ylabel('Rows (neural space dimensions)');
title('$U_m$', 'interpreter','latex');
% Third subplot: U_n_tilde
subplot(2,2,3);
imagesc(U_n_tilde);
axis image;
xlabel('Columns (basis vectors)');
ylabel('Rows (neural space dimensions)');
title('$\tilde{U}_n$', 'interpreter','latex');
% Fourth subplot: U_m_tilde
subplot(2,2,4);
imagesc(U_m_tilde);
axis image;
xlabel('Columns (basis vectors)');
ylabel('Rows (neural space dimensions)');
title('$\tilde{U}_m$', 'interpreter','latex');

% Plot hyperplanes in original neural space
figure();
us = {U_n_tilde, U_m_tilde};
for currDay = 1:2
    subplot(2,1,currDay);
    hold on;    
    % Plot data
    plot3(days{currDay}(1,:), days{currDay}(2,:), days{currDay}(3,:), 'x-');
    % Plot basis vectors
    for currDir = 1:2
        quiver3(0, 0, 0, ...
            max(abs(days{currDay}(:)))*us{currDay}(cells{currDay}(1), currDir), ...
            max(abs(days{currDay}(:)))*us{currDay}(cells{currDay}(2), currDir), ...
            max(abs(days{currDay}(:)))*us{currDay}(cells{currDay}(3), currDir), 0, ...
            'LineWidth', 4);
    end
    % Plot plane spanned by first two principal directions
    fmesh(@(s,t) us{currDay}(cells{currDay}(1),1)*s+us{currDay}(cells{currDay}(1),2)*t, ...
        @(s,t) us{currDay}(cells{currDay}(2),1)*s+us{currDay}(cells{currDay}(2),2)*t, ...
        @(s,t) us{currDay}(cells{currDay}(3),1)*s+us{currDay}(cells{currDay}(3),2)*t, ...
        [-1, 1])
    alpha(0.5);
    hold off;
    % Set plot layout properties
    legend('Data', 'Basis vector 1',...
        'Basis vector 2', 'Manifold/subpace');
    xlabel(['Neuron ' num2str(cells{currDay}(1)) ' activity']);
    ylabel(['Neuron ' num2str(cells{currDay}(2)) ' activity']);
    zlabel(['Neuron ' num2str(cells{currDay}(3)) ' activity']);
    xlim([-1,1]);
    ylim([-1,1]);
    zlim([-1,1]);
    view(75, 30);
    grid on;
    title(['Manifold in neural space of day ' num2str(currDay)]);
end

%% Section 3: principal angles between manifolds
% In this section we'll find the angles between the manifolds of the different
% recording days using another svd, now of the inner products of basis vectors

% Calculate inner products between basis vectors in U_n_tilde and U_m_tilde
innerProducts = U_n_tilde' * U_m_tilde;

% Do another svd on the inner product between basis vectors
[U_ip, S_ip, V_ip] = svd(innerProducts);

% Calculate the angles between hyperplanes from the diagonal of S_ip, which
% holds the ordered cosines of principal angles
angles = rad2deg(acos(diag(S_ip)));
disp(['Angles between hyperplanes are '...
    num2str(angles(1), 2) ', ' num2str(angles(2), 2)]);

% Sadly we can't actually plot the hyperplanes to inspect the angles:
% they live in the D-dimensional (D=4) space of the union of neurons recorded
% on day 1 and day 2. Just to illustrate the point, make two dummy basis matrices
% containing two 3D basis vectors each, and get the SVD of the inner products.

% First basis consists of x and y vector
u = [[1, 0, 0]', [0, 1, 0]'];
% Second basis matrix is first basis, but rotated 30 degrees along x
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
v = rotx(deg2rad(30)) * u;

% Do SVD, and get angles from ordered cosines on diagonal of S
[~, S, ~] = svd(u'*v);
uv = rad2deg(acos(diag(S)));

% Plot basis vectors and hyperplanes
figure()
hold on;
% Plot basis vectors u
for currDir = 1:2
    quiver3(0, 0, 0, ...
        u(1, currDir), ...
        u(2, currDir), ...
        u(3, currDir), 0, ...
        'LineWidth', 4);
end
% Plot basis vectors v
for currDir = 1:2
    quiver3(0, 0, 0, ...
        v(1, currDir), ...
        v(2, currDir), ...
        v(3, currDir), 0, ...
        'LineWidth', 4);
end
% Plot plane spanned by basis vectors u
fmesh(@(s,t) u(1,1)*s+u(1,2)*t, ...
    @(s,t) u(2,1)*s+u(2,2)*t, ...
    @(s,t) u(3,1)*s+u(3,2)*t, ...
    [-1, 1])
alpha(0.5);
% Plot plane spanned by basis vectors v
fmesh(@(s,t) v(1,1)*s+v(1,2)*t, ...
    @(s,t) v(2,1)*s+v(2,2)*t, ...
    @(s,t) v(3,1)*s+v(3,2)*t, ...
    [-1, 1])
alpha(0.5);
hold off;
% Set plot layout properties
view(60,30);
legend('Basis vector u1', 'Basis vector u2', 'Basis vector v1', 'Basis vector v2',...
    'Hyperplane spanned by u', 'Hyperplane spanned by v');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
title({'Illustration of hyperplane angle between basis sets',...
    ['Angles between hyperplanes: ' num2str(uv(1)) ', ' num2str(uv(2))]});
