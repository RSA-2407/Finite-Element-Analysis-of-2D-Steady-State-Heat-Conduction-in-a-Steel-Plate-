clc;
clear;

% Define the domain dimensions
W = 1;
H = 1;

% Number of nodes in each direction
Nx = 2;
Ny = 2;

% Element size
dx = W / (Nx-1);
dy = H / (Ny-1);

% Material properties
T_val = 25;
k_val = 50;

% Number of triangular elements
n = 2*(Nx-1)*(Ny-1);

% Allocate memory for mesh nodes
p = zeros(2, Nx*Ny);
el = zeros(3, n);

% nodes
index = 0;
for i = 1:Ny
    y = (i-1) * dy;
    for j = 1:Nx
        x = (j-1) * dx;
        index = index + 1;
        p(:, index) = [x; y];
    end
end

% elements
index = 0;
for i = 1:Ny-1
    for j = 1:Nx-1
        index1 = j + (i-1)*Nx;
        index2 = index1 + 1;
        index3 = index2 + Nx;
        index4 = index1 + Nx;

        index = index + 1;
        el(:,index) = [index1; index2; index4];

        index = index + 1;
        el(:,index) = [index2; index3; index4];
    end
end

% Plot mesh
figure;
hold on;
patch('faces', el', 'vertices', p', 'facecolor', 'c', 'edgecolor', 'k');
plot(p(1,:), p(2,:), 'o', 'MarkerFaceColor', 'r');
axis off;
hold off;

% Stiffness matrix assembly
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
for i = 1:n
    nodes = el(:,i);
    x = p(1,nodes);
    y = p(2,nodes);
    
    A0 = 0.5 * abs((x(2)-x(1))*(y(3)-y(1)) - (y(2)-y(1))*(x(3)-x(1)));
    
    A1 = 0.5 * abs((x(2)-xc)*(y(3)-yc) - (y(2)-yc)*(x(3)-xc));
    A2 = 0.5 * abs((x(1)-xc)*(y(3)-yc) - (y(1)-yc)*(x(3)-xc));
    A3 = 0.5 * abs((x(2)-xc)*(y(1)-yc) - (y(2)-yc)*(x(1)-xc));
    
    N1 = A1/A0;
    N2 = A2/A0;
    N3 = A3/A0;
    
    c = q;
    q_global(nodes(1)) = q_global(nodes(1)) + c * N1;
    q_global(nodes(2)) = q_global(nodes(2)) + c * N2;
    q_global(nodes(3)) = q_global(nodes(3)) + c * N3;
end

% Apply boundary conditions
Known = find(p(2,:) == 1); %assuming left side to be 4th side
Unknown = setdiff(1:length(p), Known);

% Partition stiffness matrix and force vector
K1 = K_global(Unknown, Unknown)
K2 = K_global(Unknown, Known)
T_calc = T_val * ones(length(Known),1);
f = q_global(Unknown) + K2 * T_calc
K2 * T_calc
q_global(Unknown)
% Solve for unknown temperatures
T_free =inv(K1)*f;

% Assemble final temperature vector
T = zeros(length(p),1);
T(Unknown) = T_free;
T(Known) = T_val;

% Display results
disp('Temperature Distribution (C):');
for i = 1:size(p,2)
    disp(['Node ', num2str(i), ' (', num2str(p(1,i)), ',', num2str(p(2,i)), '): ', num2str(T(i))]);
end

% Plot contour of temperature
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

% Plot temperature along y = 0.5
x_vals = linspace(0, W, 50);
T_vals = griddata(p(1,:), p(2,:), T, x_vals, repmat(0, size(x_vals)), 'cubic');
figure;
plot(x_vals, T_vals, '-', 'LineWidth', 2);
title('Temperature Variation along y = 0.5');
xlabel('X');
ylabel('Temperature (C)');
grid on;