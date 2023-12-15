clear all; close all; 

%% Given parameters
E = 2.1e5;
nu = 0.3;

C = E / (1 - nu^2) * [1, nu, 0;
                      nu, 1, 0;
                      0,  0, (1 - nu) / 2];

%% Model definition
Lx = 1.5;
Ly = 2/3;
nodes = [0,0; Lx/2,0; Lx,0;
         0,Ly/2; Lx/2,Ly/2; Lx*1.05,Ly/2; 
         0,Ly; Lx/2,Ly; Lx*1.2,Ly*3/4;   1.5*Lx,Ly*7/8;  2*Lx,Ly; 
         0,1.5*Ly; Lx/2,1.5*Ly; Lx,1.5*Ly; Lx*1.5,1.5*Ly; 2*Lx,1.5*Ly; 
         0,2*Ly; Lx/2,2*Ly; Lx,2*Ly; Lx*1.5,2*Ly; 2*Lx,2*Ly]; 
elements = [1,3,9,7,2,6,8,4,5;
            7,9,19,17,8,14,18,12,13;
            9,11,21,19,10,16,20,14,15];

plotElements(nodes, elements, 'Model'); grid on;
plotQuadraticElements(nodes, elements); grid on;

% Helper variables
nNodes = size(nodes, 1);
dofPerNode = 2;

nElements = size(elements, 1);
nodesPerElement = size(elements, 2);
dofPerElement = dofPerNode * nodesPerElement;
elementDOFs = globalElementDOFs(elements, 2);

%% Generalized stiffness assembly
% Set up system matrices
K = zeros(dofPerNode * nNodes);
u = zeros(dofPerNode * nNodes, 1);
f = zeros(dofPerNode * nNodes, 1);

% Set up stiffness matrix by looping over elements
for e = 1 : nElements
    KE = elementStiffness(C, elementNodes, integrationOrder);
    
    % Use element connectivity to add element stiffness to global stiffness matrix
end

% Load
f = f;

% Solve system
fixeddofs = [];
freedofs = setdiff(1 : nNodes * dofPerNode, fixeddofs);

u(freedofs) = K(freedofs, freedofs) \ f(freedofs);

% Plot displacement as x0 + u
scale = 3;
uNodes = transpose(reshape(u, dofPerNode, nNodes));
uMagnitude = vecnorm(uNodes, 2, 2);
% plotElements(nodes + scale * uNodes, elements, '|u|', uMagnitude);
plotQuadraticElements(nodes + scale * uNodes, elements);

%% Additional functions
function elementDOFs = globalElementDOFs(elements, dim)
    % Return global DOFs for elements
    % INPUT
    %   elements (array)    nodes per element as (nElements x nodesPerElement) array
    %   dim (scalar)        dimension
    % OUTPUT
    %   elementDOFs (array) global DOFs per element as (nElements x (dim * nodesPerElement)) array

    switch dim
        case 1
            elementDOFs = elements;
        case 2
            nodes = reshape(elements', numel(elements), 1);
            elementDOFs = [2 * nodes - 1, 2 * nodes]';
            elementDOFs = transpose(reshape(elementDOFs(:), 2 * size(elements, 2), size(elements, 1)));
        otherwise
            error("not implemented");
    end
end

function KE = elementStiffness(C, nodes, integrationOrder)
    % Compute element stiffness matrix
    % INPUT
    %   C (array)                   Constitutive matrix
    %   nodes (array)               Nodes for element
    %   integrationOrder (scalar)   Integration order for gauss quadrature
    % OUTPUT
    %   KE (array)                  Element stiffness matrix

end

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

function [h, dh] = quadraticShapeFunctions(r, s)
    % Evaluate second order shape functions for quadrilateral element
    % INPUT
    %   r (array)   r coordinate for evaluation of shape functions and derivatives in natural coordinates
    %   s (array)   s coordinate for evaluation of shape functions and derivatives in natural coordinates
    % OUTPUT
    %   h (array)   Shape functions at (r, s), row per given point
    %   dh (array)  Derivatives of shape functions at first (r, s)

    h = [];
    dh = [];
end

% H and B matrices
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

function B = bMatrix(dhdx, order)
    % Evaluate B matrix for quadrilaterial element
    % INPUT
    %   r (array)   r coordinate for for evaluation of shape functions and derivatives in natural coordinates
    %   s (array)   s coordinate for for evaluation of shape functions and derivatives in natural coordinates
    % OUTPUT
    %   B (array)   B matrix for linear quadrilateral element in natural coordinates

    % Initialize B matrix
    if order == 1
        nNodes = 4;

    else
        nNodes = 9;
    end
    B = zeros(3, 2 * nNodes);

    % Define linear indices for first B matrix -> B1 corresponding to node 1
    indicesBi = [1; 5; 3; 6];
    indicesdHi = [1; 2; 2; 1];

    % Repeat indices for node 2 to node X
    indicesB = indicesBi + (0 : 6 : 6 * nNodes - 6);
    indicesdH = indicesdHi + (0 : 2 : nNodes * 2 - 2);

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
        % Reorder nodes for second order element for Matlab patch
        elementNodes = elementNodes([1 5 2 6 3 7 4 8], :);
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

function plotQuadraticElements(nodes, elements)
    % Plot second-order elements by computing points on edge of element in (x, y)
    % INPUT
    %   nodes (array)       Nodal coordinates as (nNodes x 2) array
    %   elements (array)    Nodal indices per element as (nElements x 4) array

    newFigure();
    hold on;

    % Set up vectors for (r,s)-coordinates along edge of isoparametric element
    rs = [-1 : 0.1 : 1]';
    r = [rs; ones(numel(rs) - 2, 1); flip(rs); -1 * ones(numel(rs) - 2, 1); -1];
    s = [-1 * ones(numel(rs), 1); rs(2 : end - 1); ones(numel(rs), 1); flip(rs(2 : end - 1)); -1];
    
    % Loop over elements to plot contour for each element
    for e = 1 : size(elements, 1)
        elementNodes = nodes(elements(e, :), :);

        % Get shape functions for points along edge
        h = quadraticShapeFunctions(r, s);

        % Compute (x, y) data for (r, s) along edges 
        x = h * elementNodes(:, 1);
        y = h * elementNodes(:, 2);

        % Connect data points by plotting lines
        plot(x, y, 'k');

        % Optional: Plot nodes
        scatter(elementNodes(:, 1), elementNodes(:, 2), 'o', 'filled', ...
            'MarkerEdgeColor','k','MarkerFaceColor','w','SizeData', 100);
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