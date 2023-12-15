clear all; close all; 

%% Given parameters
E = 2.1e5;
nu = 0.3;

C = E / (1 - nu^2) * [1, nu, 0;
                      nu, 1, 0;
                      0,  0, (1 - nu) / 2];

%% Model definition
nodes = [-1, -1; 0.8, -1.2; 0.6, 1; -1.2, 0.6];
elements = [1, 2, 3, 4];

plotElements(nodes, elements, 'Model'); grid on;

% Helper variables
nNodes = size(nodes, 1);
dofPerNode = 2;

nElements = size(elements, 1);
nodesPerElement = size(elements, 2);
dofPerElement = dofPerNode * nodesPerElement;

% Look at integration points in global coordinate system
[gauss_rs, wg] = gaussianPoints2D(2);
gauss_xy = zeros(size(gauss_rs));
for g = 1 : size(gauss_rs, 1)
    % Evaluate shape function derivatives at gauss point
    [h, ~] = linearShapeFunctions(gauss_rs(g, 1), gauss_rs(g, 2));
    
    % for j = 1 : nodesPerElement
    %     gauss_xy(g, :) = gauss_xy(g, :) + h(j) * nodes(j, :);
    % end
    % Shorter:
    gauss_xy(g, :) = h * nodes;
end
hold on; scatter(gauss_xy(:, 1), gauss_xy(:, 2), Marker = '+', SizeData = 120, LineWidth = 2); hold off;

% Compute center
center = linearShapeFunctions(0, 0) * nodes;
hold on; scatter(center(1), center(2), Marker = '+', SizeData = 120, LineWidth = 2, color = 'r'); hold off;

% Stiffness matrix
Ke = zeros(dofPerNode * nodesPerElement);
for g = 1 : size(gauss_rs, 1)
    [~, dh] = linearShapeFunctions(gauss_rs(g, 1), gauss_rs(g, 2));

    % Jacobian at gauss point
    J = dh * nodes;
    detJ = det(J);
    invJ = inv(J);

    % % Computation of submatrices for stiffness matrix
    % for i = 1 : size(nodes, 1)
    %     dhidx = invJ * dh(:, i);
    %     % dhidx = J \ dh(:, 1);
    %     Bi = [dhidx(1), 0;
    %           0, dhidx(2);
    %           dhidx(2), dhidx(1)];
    %     for j = i : size(nodes, 1)
    %         dhjdx = invJ * dh(:, j);
    %         Bj = [dhjdx(1), 0;
    %             0, dhjdx(2);
    %             dhjdx(2), dhjdx(1)];
    %         Kij = Bj' * C * Bi * abs(detJ);
    %     end
    % end
    % Shorter:
    dhdx = invJ * dh;
    B = bMatrix(dhdx);
    Ke = Ke + B' * C * B * abs(detJ);
end

% Body load
f0 = [];
fbe = [];
for g = 1 : size(gauss_rs, 1)
    % Evaluate shape function derivatives at gauss point
    [h, dh] = linearShapeFunctions();

    % Set up jacobian matrix
    J = [];

    % Evaluation of summand for body load integration
    for j = 1 : nodesPerElement
        fbe = fbe;
    end
end

% Edge load
fse = [];

%% Additional functions
% Shape functions
function [h, dh] = linearShapeFunctions(r, s)
    % Evaluate shape functions for linear quadrilateral element
    % INPUT
    %   r (array)   r coordinate for evaluation of shape functions and derivatives in natural coordinates
    %   s (array)   s coordinate for evaluation of shape functions and derivatives in natural coordinates
    % OUTPUT
    %   h (array)           Shape functions at (r, s), row per given point
    %   dh (array)          Derivatives of shape functions at first (r, s)
    h1 = 1/4 * (1 - r) .* (1 - s);
    h2 = 1/4 * (1 + r) .* (1 - s);
    h3 = 1/4 * (1 + r) .* (1 + s);
    h4 = 1/4 * (1 - r) .* (1 + s);

    dh1dr = -1/4 * (1 - s);
    dh1ds = -1/4 * (1 - r);
    dh2dr = 1/4 * (1 - s);
    dh2ds = -1/4 * (1 + r);
    dh3dr = 1/4 * (1 + s);
    dh3ds = 1/4 * (1 + r);
    dh4dr = -1/4 * (1 + s);
    dh4ds = 1/4 * (1 - r);

    h = [h1, h2, h3, h4];
    dh = [dh1dr, dh1ds;
          dh2dr, dh2ds;
          dh3dr, dh3ds;
          dh4dr, dh4ds]';
end

% H and B matrices and Jacobian
function H = hMatrix(r, s)
    % Evaluate H matrix for linear quadrilaterial element
    % INPUT
    %   r (array)   r coordinate for for evaluation of shape functions and derivatives in natural coordinates
    %   s (array)   s coordinate for for evaluation of shape functions and derivatives in natural coordinates
    % OUTPUT
    %   H (array)   B matrix for linear quadrilateral element in natural coordinates

    h = linearShapeFunctions(r, s);
    H = [h(1), 0, h(2), 0, h(3), 0, h(4), 0;
         0, h(1), 0, h(2), 0, h(3), 0, h(4)];
end

function B = bMatrix(dhdx)
    % Evaluate B matrix for linear quadrilaterial element
    % INPUT
    %   r (array)   r coordinate for evaluation of shape functions and derivatives in natural coordinates
    %   s (array)   s coordinate for evaluation of shape functions and derivatives in natural coordinates
    % OUTPUT
    %   B (array)   B matrix for linear quadrilateral element in natural coordinates

    % Initialize B matrix
    B = zeros(3, 8);

    % Define linear indices for first B matrix -> B1 corresponding to node 1
    indicesBi = [1; 5; 3; 6];
    indicesdHi = [1; 2; 2; 1];

    % Repeat indices for node 2 to node 4
    indicesB = indicesBi + (0 : 6 : 6 * 4 - 6);
    indicesdH = indicesdHi + (0 : 2 : 4 * 2 - 2);

    % Assign dh entries to B matrix
    B(indicesB) = dhdx(indicesdH);
end

% Plot functions
function p = plotElements(nodes, elements, title, color)
    % Plot elements by plotting a patch for each element
    % INPUT
    %   nodes (array)       Nodal coordinates as (nNodes x 2) array
    %   elements (array)    Nodal indices per element as (nElements x 4) array
    %   color (array)       OPTIONAL: Color for nodes or elements per row
    % OUTPUT
    %   p (array)           Patches for elements

    newFigure(title);
    hold on;
    for e = 1 : size(elements, 1)
        elementNodes = nodes(elements(e, :), :);
        p(e) = patch(elementNodes(:, 1), elementNodes(:, 2), 'w');
        if nargin > 3
            if size(color, 1) == size(nodes, 1)
                % Interpolate between nodal values
                p(e).set('FaceColor', 'interp', 'FaceVertexCData', color(elements(e, :)));
            elseif size(color, 1) == size(elements, 1)
                % Color per element
                p(e).set('FaceColor', 'flat', 'FaceVertexCData', color(e));
            end
            clim([min(color), max(color)]);
            % p(e).EdgeAlpha = 0;
        else
            p(e).FaceAlpha = 0;
        end
    end
    hold off;
end

function f = newFigure(figureTitle)
    % Open figure and set several common options
    % INPUT
    %   figureTitle (string/char array)     OPTIONAL: title for figure
    % OUTPUT
    %   f (figure) 

    f = figure;
    f.Position([3, 4]) = 1.5 * f.Position([3, 4]);
    f.Position(2) = f.Position(2) - f.Position(4)/3;
    axis equal; 
    colormap jet;
    colorbar;
    if nargin > 0
        title(figureTitle);
    end
end