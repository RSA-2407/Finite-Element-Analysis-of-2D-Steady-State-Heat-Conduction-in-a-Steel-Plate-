clc;
clear;

% Define the domain dimensions
W = 1;
H = 1;

% Number of nodes in each direction
Nx = 6;
Ny = 6;

% Element size
dx = W / (Nx-1);
dy = H / (Ny-1);

% Material properties
T_val = 25;
k_val = 50;

% nodes
p = zeros(2, Nx*Ny);
index = 0;
for i = 1:Ny
    for j = 1:Nx
        index = index + 1;
        p(:, index) = [(j-1)*dx; (i-1)*dy];
    end
end

% Center nodes
x_cells = Nx-1;
y_cells = Ny-1;
centers = zeros(2, x_cells*y_cells);
index = 0;
for i = 1:y_cells
    for j = 1:x_cells
        index = index + 1;
        centers(:, index) = [(j-0.5)*dx; (i-0.5)*dy];
    end
end
p = [p centers]; % Combine original and center nodes

% Generate triangular elements
el = [];
center_offset = Nx * Ny;
for i = 1:y_cells
    for j = 1:x_cells
        % Corner node indices
        index1 = (i-1)*Nx + j;
        index2 = index1 + 1;
        index3 = index2 + Nx;
        index4 = index3 - 1;
        
        % Center node index
        center_idx = center_offset + (i-1)*x_cells + j;
        
        % Create four triangular elements per cell ensuring symmetry
        el = [el [index1; index2; center_idx]];
        el = [el [index2; index3; center_idx]];
        el = [el [index3; index4; center_idx]];
        el = [el [index4; index1; center_idx]];
    end
end

figure;
hold on;
patch('faces', el', 'vertices', p', 'facecolor', 'c', 'edgecolor', 'k');
plot(p(1,:), p(2,:), 'o', 'MarkerFaceColor', 'r');
axis off;
hold off;

% Stiffness matrix assembly
n = length(el);
K_global = zeros(length(p));

for i = 1:n
    nodes = el(:,i);
    x = p(1, nodes);
    y = p(2, nodes);
    
    J = [x(1)-x(3), x(2)-x(3); y(1)-y(3), y(2)-y(3)];
    A = 0.5 * abs(det(J));
    Q = [1, 0, -1; 0, 1, -1];
    K = k_val * A * (Q' / J) * (Q' / J)';
    
    for j = 1:3
        for k = 1:3
            K_global(nodes(j), nodes(k)) = K_global(nodes(j), nodes(k)) + K(j,k);
        end
    end
end

% Load vector (heat source)
xc = 0.5;
yc = 0.5;
q = 1000;
q_global = zeros(size(p,2),1);

% Distribute heat source among elements
for i = 1:length(p)
    x = p(1,i);
    y = p(2,i);
    if x == 0.5 && y == 0.5
        q_global(i) = 1000;
        break;
    end
end

% Apply boundary conditions
Known = find(p(1,:) == 0); %assuming left side to be 4th side
Unknown = setdiff(1:length(p), Known);

% Partition stiffness matrix and force vector
K1 = K_global(Unknown, Unknown);
K2 = K_global(Unknown, Known);
T_calc = T_val * ones(length(Known),1);
f = q_global(Unknown) - K2 * T_calc;

% Solve for unknown temperatures
T_free = K1 \ f;

% Assemble final temperature vector
T = zeros(length(p),1);
T(Unknown) = T_free;
T(Known) = T_val;

% Display results
disp('Temperature Distribution (C):');
for i = 1:size(p,2)
    disp(['Node ', num2str(i), ' (', num2str(p(1,i)), ',', num2str(p(2,i)), '): ', num2str(T(i))]);
end

% Plot temperature contour
grid_x = linspace(0, W, 50);
grid_y = linspace(0, H, 50);
[Xq, Yq] = meshgrid(grid_x, grid_y);
Tq = griddata(p(1,:), p(2,:), T, Xq, Yq, 'cubic');
figure;
contourf(Xq, Yq, Tq, 200, 'LineColor', 'none');
colorbar;
title('Temperature Contour Plot');
xlabel('X');
ylabel('Y');

% Plot temperature along horizontal line
x_vals = linspace(0, W, 50);
T_vals = griddata(p(1,:), p(2,:), T, x_vals, repmat(0.5, size(x_vals)), 'cubic');
figure;
plot(x_vals, T_vals, '-', 'LineWidth', 2);
title('Temperature Variation along y = 0.5');
xlabel('X');
ylabel('Temperature (C)');
grid on;