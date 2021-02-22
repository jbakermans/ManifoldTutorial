% PART 1 - https://www.crowdcast.io/e/sara-sollas-world-wide, 15:00-25:00
% Load synthetic data and find low-dimensional manifold.
% Project data on neural modes; see how latent trajectory generates data.

% This script has a number of blanks like this: '___' (triple underscore).
% You can try filling them out and then running the result section by section.
% Press cmd + Enter to have Matlab run the current section
% Solutions can be found in part1_solutions.m

% Jacob Bakermans, February 2021
close all;
clear all;
%% Section 1: Load data
% In this section we'll generate some synthetic neural data for 3 neurons,
% and plot the data over time and as a trajectory in neuron activity space

% Generate synthetic data
X_D = generate_data_1(); % Data matrix: N_neurons x N_timebins

% Plot data: firing rate for each neuron
figure(); 
% First subplot: neuron firing rate through time
subplot(1,2,1);
hold on;
for currNeuron = 1:size(X_D,1)
    plot(X_D(currNeuron, :));
end
hold off;
% Set plot layout properties
legend('Neuron 1', 'Neuron 2', 'Neuron 3');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Neuron activity through time');
% Second subplot: trajectory in neural space
subplot(1,2,2);
% Plot data
plot3(X_D(1,:), X_D(2,:), X_D(3,:), 'x-');
% Set plot layout properties
xlabel('Neuron 1 activity');
ylabel('Neuron 2 activity');
zlabel('Neuron 3 activity');
view(75, 30);
grid on;
title('Trajectory in neural space');

%% Section 2: Do PCA to get neural modes/principal components 
% In this section we'll find the flat subspace/linear manifold in the data by PCA.
% This manifold is spanned by the first two principal components, or neural modes.

% A PCA is an eigendecomposition of the covariance matrix of a series of datapoints
% Remove mean firing rate for every neuron
X_D = '___';
% Calculate covariance matrix between neurons
cov_dat = '___';
% Do eigenvalue decomposition of covariance
[V, D] = '___';

% V has eigenvectors in columns: these are the principal components, or neural modes
% D has corresponding eigenvalues on diagonal
% The higher the eigenvalue, the more variance of activity is explained by that neural mode
% Therefore it's useful to sort modes by eigenvalue: highest eigenvalues first

% Sort eigenvectors (columns) in V by descending eigenvalue (diagonal) in D
V = '___';

% Plot neural trajectory and principal components
figure();
hold on;
% Plot data
plot3(X_D(1,:), X_D(2,:), X_D(3,:), 'x-');
% Plot principal components
for currDir = 1:3
    quiver3(0, 0, 0, max(abs(X_D(:)))*V(1, currDir), ...
        max(abs(X_D(:)))*V(2, currDir), max(abs(X_D(:)))*V(3, currDir), 0, ...
        'LineWidth', 4);
end
% Plot plane spanned by first two principal components
fmesh(@(s,t)V(1,1)*s+V(1,2)*t, @(s,t)V(2,1)*s+V(2,2)*t, @(s,t)V(3,1)*s+V(3,2)*t, ...
    [-1, 1])
alpha(0.5);
hold off;
% Set plot layout properties
legend('Data', 'Principal component 1',...
    'Principal component 2', 'Principal component 3', 'Manifold/subpace');
xlabel('Neuron 1 activity');
ylabel('Neuron 2 activity');
zlabel('Neuron 3 activity');
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);
view(75, 30);
grid on;

%% Section 3: Go from neuron activity space to neural mode space
% In this section, instead of plotting the trajectory in neuron activity space, we'll move
% the trajectory to the neural mode space: the space spanned by the neural modes
% (u and v on the slides). We'll get what Solla calls latent rates. You could think about
% latent rates as generating the measured neuron activity by linear combination.

% The first two neural modes span the manifold. Select the corresponding columns from V
V_tilde = '___';

% Now move trajectory in neuron activity space to neural mode space to get latent rates
L = '___';

% Plot trajectory in neural mode space
figure(); 
plot(L(1,:), L(2,:), '-x');
axis square;
xlabel('Neural mode 1 (u)');
ylabel('Neural mode 2 (v)');
title('Latent variables: data projected on manifold, in neural mode space');

% Solla's "alternative interpretation": latent rates generate neural activity
% Use linear combinations of the latent trajectory to recover original activity
X_gen = '___';
% Plot neuron activity through time, and compare to latent activity and recovered activity
figure(); 
% First subplot: original data
subplot(3,1,1)
hold on;
for currNeuron = 1:size(X_D,1)
    plot(X_D(currNeuron, :));
end
hold off;
% Set plot layout properties
legend('Neuron 1', 'Neuron 2', 'Neuron 3');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Neuron activity through time');
% Second subplot: latent activity
subplot(3,1,2)
hold on;
for currNeuron = 1:size(L,1)
    plot(L(currNeuron, :));
end
hold off;
% Set plot layout properties
legend('Latent unit 1', 'Latent unit 2');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Latent activity through time');
% Third subplot: recovered original data from latent units
subplot(3,1,3)
hold on;
for currNeuron = 1:size(X_gen,1)
    plot(X_gen(currNeuron, :));
end
hold off;
% Set plot layout properties
legend('Neuron 1', 'Neuron 2', 'Neuron 3');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Recovered activity through time');

%% Section 5: peek ahead: SVD
% In this section we'll get a head start with the next bit of the talk,
% where the Singular Value Decomposition (SVD) is going to show up many times.
% Turns out that's basically the same as what we've been doing so far...

% Do a singular value decomposition of the input data X_D, 
% so that X_D = U_svd * S_svd * V_svd'
[U_svd, S_svd, V_svd] = '___';

% U_svd is a matrix with the eigenvectors of X*X' in its columns;
% V_svd is a matrix with the eigenvectors of X'*X in its columns.
% That means that U_svd should be the same as our V in the previous sections,
% except for sign flips: the same eigenvectors may point in the opposite direction

% Plot them next to each other to see if that's true
figure()
% First subplot: U from SVD
subplot(1,3,1);
imagesc(U_svd);
axis square;
title('U from SVD');
% Second subplot: V from PCA
subplot(1,3,2);
imagesc(V);
axis square;
title('V from PCA');
% Third subplot: inner products between U and V. 
% All values on the diagonal should be 1 (parallel) or -1 (opposite).
subplot(1,3,3);
imagesc(U_svd'*V);
axis square;
title('Column inner products');
colorbar;