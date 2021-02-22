% PART 3 - https://www.crowdcast.io/e/sara-sollas-world-wide, 55:00-1:00:00
% Use basis vectors for different days in common ambient space to align
% latent dynamics via Canonical Correlation Analysis, using QR-decomposition
% Jacob Bakermans, February 2021
close all;
clear all;
%% Section 1: Prepare data (copy from part 2)
% In this section, we'll repeat the steps from part 2: load data, combine
% into space of union of neurons across days, do svd, and find basis for 
% low-dimensional manifold for each day

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

% Do singular value decomposition for both days. U contains basis vectors in
% D dimensions in its columns and provides a basis in the space of neurons.
[U_n, S_n, V_n] = svd(X_n);
[U_m, S_m, V_m] = svd(X_m);

% Set 'flat dimension' d
d = 2;

% Select first d basis vectors from U_n and U_m
U_n_tilde = U_n(:, 1:d);
U_m_tilde = U_m(:, 1:d);

%% Section 2: project data on manifolds to get latent variables
% In this section, we'll project the data of both days on the low-dimensional
% manifold for the corresponding day, to get trajectories in space of latent
% variables/neural modes. These are not yet aligned.

% Project data on d-dimensional manifolds to get latent variables
L_n = U_n_tilde' * X_n;
L_m = U_m_tilde' * X_m;

% Calculate unaligned canonical correlations from inner products
CC_unaligned = L_n*L_m';

% Plot resulting latent trajectories
% Plot hyperplanes in original neural space
figure();
us = {U_n_tilde, U_m_tilde};
ls = {L_n, L_m};
for currDay = 1:2
    % First subplot for each day: trajectory in original neural activity space
    subplot(2,2,(currDay-1)*2+1);    
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
    % Second subplot: latent variables/trajectory in neural modes space
    subplot(2,2,(currDay-1)*2+2);
    plot(ls{currDay}(1,:), ls{currDay}(2,:), '-x');
    % Set plot layout properties
    xlabel('Neural mode 1');
    ylabel('Neural mode 2');
    title(['Latent variables on day ' num2str(currDay)]);
end

%% Section 3: Canonical Correlation Analysis (CCA)
% In this section we'll do canonical correlation analysis by QR decomposition
% of the latent variables, and svd of the orthonormal latent variables.

% If you transpose the latent variables, you'll get a T x d data matrix, 
% with the latent activity at every tim bin in the rows. 
% The QR decomposition of L' finds an orthonormal basis for L', 
% where Q holds d T-dimensional orthonormal vectors and R the coefficients
% of the linear combination that will recover input L'.

% Do the QR decomposition of transposed latent variables. Hint: to get T x d matrix 
% of orthornormal basis vectors, do "economy size" QR decomposition
[Q_n, R_n] = qr(L_n', 0);
[Q_m, R_m] = qr(L_m', 0);

% Do SVD of inner product of Qs to get aligned canonical correlation coefficients
[U_q, S_q, V_q] = svd(Q_n' * Q_m);

% Q being orthornomal means that the inner product of columns i and j of Q 
% is 1 if i == j, and 0 if i ~= j. Also L' should be recoverd by Q*R.
% Make a plot to verify both are true.
figure()
qs = {Q_n, Q_m};
rs = {R_n, R_m};
for currDay = 1:2
    % First subplot for each day: latent variables  
    subplot(2,3,(currDay-1)*3+1);    
    plot(ls{currDay}(1,:), ls{currDay}(2,:), '-x');
    % Set plot layout properties
    xlabel('Neural mode 1');
    ylabel('Neural mode 2');
    title(['Latent variables on day ' num2str(currDay)]);
    % Second subplot for each day: QR decomposition of latent variables
    % Expect inner products to be identity matrix, and different to be 0
    subplot(2,3,(currDay-1)*3+2);    
    imagesc(Q_n' * Q_n);
    % Set plot layout properties    
    xlabel('Basis vectors');
    ylabel('Basis vectors');
    title({'Inner product between columns of Q',...
        ['Maximum of L_n^T - Q_n*R_n: ' ...
        num2str(max(max(abs((Q_n * R_n) - L_n'))),2)]});
    % Second subplot for each day: latent variables in orthonormal basis 
    subplot(2,3,(currDay-1)*3+3);    
    plot(qs{currDay}(:,1), qs{currDay}(:,2), '-x');
    % Set plot layout properties
    xlabel('Orthonormal mode 1');
    ylabel('Orthonormal mode 2');
    title(['Orthonormal latent variables on day ' num2str(currDay)]);  
end

%% Section 4: Align latent variables
% In this section we'll do the final step in aligning latent dynamics: we'll
% use the outcomes of CCA to calculate and apply the linear transformation
% that takes the latent variables in the neural modes of different days into
% a single aligned manifold.

% Calculate transformation matrices M that move latent variables to manifold 
% directions in which latent variables are maximally correlated across days
M_n = inv(R_n) * U_q;
M_m = inv(R_m) * V_q;

% Use transformation matrices M to carry out the linear transformation of 
% latent variables into the aligned space
L_n_aligned = M_n' * L_n;
L_m_aligned = M_m' * L_m;

% Finally: plot results of aligned latent variables in space of maximally
% correlated neural modes found via canonical correlation analysis
figure()
for currDay = 1:2
    % First subplot for each day: latent variables  
    subplot(2,4,(currDay-1)*4+1);    
    plot(ls{currDay}(1,:), ls{currDay}(2,:), '-x');
    % Set plot layout properties
    xlabel('Neural mode 1');
    ylabel('Neural mode 2');
    title(['Latent variables on day ' num2str(currDay)]);
    % Second subplot for each day: latent variables in orthonormal basis 
    subplot(2,4,(currDay-1)*4+2);    
    plot(qs{currDay}(:,1), qs{currDay}(:,2), '-x');
    % Set plot layout properties
    xlabel('Orthonormal mode 1');
    ylabel('Orthonormal mode 2');
    title(['Orthonormal latent variables on day ' num2str(currDay)]);  
end
% One big combined subplot: both latent variables in the same aligned space
subplot(2,4,[3,4,7,8]);
hold on;
plot(L_n_aligned(1,:), L_n_aligned(2,:), '-x');
plot(L_m_aligned(1,:), L_m_aligned(2,:), '-x');
hold off;
% Set plot layout properties
legend('Day 1', 'Day 2');
xlabel('Aligned mode 1');
ylabel('Aligned mode 2');
title('Aligned latent variables across days');

%% Section 5: Post scriptum - was QR decomposition necessary?
% In this section we'll arrive at the same result taking an alternative
% route: directly use the leading factors from the initial SVD of the data,
% instead of calculating the QR decomposition of the latent variables. The 
% factors (columns of V in the initial SVD; principal components across time)
% are orthonormal already. What is the advantage of doing QR decomposition?

% Select first d basis vectors from V_n and V_m. These are 'factors', 
% vectors of linear combinations across time bins instead of neurons.
V_n_tilde = V_n(:, 1:d);
V_m_tilde = V_m(:, 1:d);

% Do SVD of inner products between factors
[U_v, S_v, V_v] = svd(V_n_tilde' * V_m_tilde);

% Now the matrix of singular values of the data is going to play to role of
% R in the QR decomposition: it holds coefficients for linear combination of
% basis vectors. Only select the first d rows and columns.
S_n_tilde = S_n(1:d, 1:d);
S_m_tilde = S_m(1:d, 1:d);

% Calculate matrices that transform data to aligned neural mode space -
% equivalent to the matrices M_m and M_n before - from SVD results
M_n_svd = inv(S_n_tilde) * U_v;
M_m_svd = inv(S_m_tilde) * V_v;

% Apply transformation to get latent variables in aligned space
L_n_aligned_svd = M_n_svd' * L_n;
L_m_aligned_svd = M_m_svd' * L_m;

% Make a plot to compare the transformation matrices from QR decomposition
% to those obtained directly from factors of original SVD
figure();
% First subplot: M_n
subplot(2,2,1);
imagesc(M_n);
axis image;
xlabel('Columns');
ylabel('Rows');
title('M_n from QR decomposition');
% Second subplot: M_m
subplot(2,2,2);
imagesc(M_m);
axis image;
xlabel('Columns');
ylabel('Rows');
title('M_m from QR decomposition');
% Third subplot: M_n_svd
subplot(2,2,3);
imagesc(M_n_svd);
axis image;
xlabel('Columns');
ylabel('Rows');
title('M_n directly from SVD');
% Fourth subplot: U_m_tilde
subplot(2,2,4);
imagesc(M_m_svd);
axis image;
xlabel('Columns');
ylabel('Rows');
title('M_m directly from SVD');

% Now plot basis vectors for alignment and aligned latent variables for
% both days, on different rows for different approaches, to compare them
figure()
vs = {V_n_tilde, V_m_tilde};
for currDay = 1:2
    % First row subplots: orthonormal latent variables for QR decomposition
    subplot(2,3,currDay);    
    plot(qs{currDay}(:,1), qs{currDay}(:,2), '-x');
    % Set plot layout properties
    xlabel('Orthonormal mode 1');
    ylabel('Orthonormal mode 2');
    title(['Orthonormal latent variables Q on day ' num2str(currDay)]); 
    % Second row subplots: factors directly from initial SVD
    subplot(2,3,3 + currDay);    
    plot(vs{currDay}(:,1), vs{currDay}(:,2), '-x');
    % Set plot layout properties
    xlabel('Factor 1');
    ylabel('Factor 2');
    title(['Factors (PCs across time) of V on day ' num2str(currDay)]);  
end
% Fourth subplot for each day: aligned latent variables using QR decomposition
subplot(2,3,3);
hold on;
plot(L_n_aligned(1,:), L_n_aligned(2,:), '-x');
plot(L_m_aligned(1,:), L_m_aligned(2,:), '-x');
hold off;
% Set plot layout properties
legend('Day 1', 'Day 2');
xlabel('Aligned mode 1');
ylabel('Aligned mode 2');
title('Aligned latent from QR decomposition');
% Fourth subplot for each day: aligned latent variables using initial SVD factors
subplot(2,3,6);
hold on;
plot(L_n_aligned_svd(1,:), L_n_aligned_svd(2,:), '-x');
plot(L_m_aligned_svd(1,:), L_m_aligned_svd(2,:), '-x');
hold off;
% Set plot layout properties
legend('Day 1', 'Day 2');
xlabel('Aligned mode 1');
ylabel('Aligned mode 2');
title('Aligned latent from initial SVD factors');
