clc;
clear;

% Define the domain dimensions
W = 1;
H = 1;

% Number of nodes in each direction
Nx = 26;
Ny = 26;

% Element size
dx = W / (Nx-1);
dy = H / (Ny-1);

% Material properties
T_val = 25;
k_val = 50;
h_val = 10;

% generating nodes
p = [];

% bottom boundary nodes
for j = 1:Nx
    p = [p [W*(j-1)/(Nx-1); 0]];
end

% right boundary nodes
for i = 2:Ny
    p = [p [W; H*(i-1)/(Ny-1)]];
end

% top boundary nodes
for j = Nx-1:-1:1
    p = [p [W*(j-1)/(Nx-1); H]];
end

% left boundary nodes
for i = Ny-1:-1:2
    p = [p [0; H*(i-1)/(Ny-1)]];
end

% center node
center_point = [W/2; H/2];
p = [p center_point];
center_idx = size(p, 2);

% Generate triangular elements
el = [];
num_boundary = center_idx - 1;

% Connect adjacent boundary nodes with the center
for i = 1:num_boundary
    next_i = i + 1;
    if next_i > num_boundary
        next_i = 1;
    end
    el = [el [i; next_i; center_idx]];
end

% Plot the mesh
figure;
hold on;
patch('faces', el', 'vertices', p', 'facecolor', 'c', 'edgecolor', 'k');
plot(p(1,1:end), p(2,1:end), 'o', 'MarkerFaceColor', 'r');
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

% Extract the left edge of the mesh (nodes with x = 0)
left_edge = find(p(1,:) == 0);

[~, s_o] = sort(p(2, left_edge));
left_edge = left_edge(s_o);

% left edge elements
l_lines = zeros(2, length(left_edge)-1);
for i = 1:length(left_edge)-1
    l_lines(:,i) = [left_edge(i); left_edge(i+1)];
end

% Add convection contribution to the global stiffness matrix
for i = 1:size(l_lines, 2)
    node1 = l_lines(1,i);
    node2 = l_lines(2,i);
    x1 = p(1, node1); y1 = p(2, node1);
    x2 = p(1, node2); y2 = p(2, node2);
    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    
    S = L/6 * [2 1; 1 2];

    K_global(node1, node1) = K_global(node1, node1) + h_val * S(1,1);
    K_global(node1, node2) = K_global(node1, node2) + h_val * S(1,2);
    K_global(node2, node1) = K_global(node2, node1) + h_val * S(2,1);
    K_global(node2, node2) = K_global(node2, node2) + h_val * S(2,2);
end 

% Load vector (heat source)
xc = 0.5;
yc = 0.5; 
q = 1000;
q_global = zeros(size(p,2),1);

% Heat source at center
for i = 1:length(p)
    x = p(1,i);
    y = p(2,i);
    if x == 0.5 && y == 0.5
        q_global(i) = 1000;
        break;
    end
end

% Flux
F_conv = zeros(size(p,2), 1);
for i = 1:size(l_lines, 2)
    node1 = l_lines(1,i);
    node2 = l_lines(2,i);
    x1 = p(1, node1); y1 = p(2, node1);
    x2 = p(1, node2); y2 = p(2, node2);
    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    
    F_vec = h_val * T_val * L/ 2 * [1; 1];
    
    F_conv(node1) = F_conv(node1) + F_vec(1);
    F_conv(node2) = F_conv(node2) + F_vec(2);
end

% Solve 
q_global = q_global + F_conv;
T = K_global \ q_global;

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