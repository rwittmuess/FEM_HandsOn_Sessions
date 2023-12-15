close all; 

% MODEL
%  6   --   7   --   8
%  |        |        |
%  |   2    |   3    |
%  3   --   4   --   5
%  |        |
%  |   1    | 
%  1   --   2
%
% ELEMENT
%  4   --   3 
%  |        |
%  |        | 
%  1   --   2

% Given parameters
lx = 1.5;
ly = 2/3;

% Material
E = 2.1e5;
nu = 0.3;
C = E / (1 - nu^2) * [1, nu, 0;
                      nu, 1, 0;
                      0, 0, (1 - nu) / 2];

% Set up nodes and elements
nodes = [0, 0;
         lx, 0;
         0, ly;
         lx, ly;
         2*lx, ly;
         0, 2*ly;
         lx, 2*ly;
         2*lx, 2*ly];
nNodes = size(nodes, 1);
dofPerNode = 2;

% Nodal indices per element
elements = [1 2 4 3;
            3 4 7 6;
            4 5 8 7];
nElements = size(elements, 1);
nodesPerElement = size(elements, 2);
dofPerElement = dofPerNode * nodesPerElement;

% Optional: Plot structure
p = plotElements(nodes, elements, '');

% Set up system matrices
K = zeros(dofPerNode * nNodes);
u = zeros(dofPerNode * nNodes, 1);
f = zeros(dofPerNode * nNodes, 1);

% Integrate stiffness by looping over elements
[xg, wg] = gaussianPoints2D();
for e = 1 : nElements
    KE = zeros(dofPerElement);
    fE = zeros(dofPerElement, 1);

    % Get element corner coordinates
    xmin = nodes(elements(e, 1),1);
    xmax = nodes(elements(e, 2),1);
    ymin = nodes(elements(e, 1),2);
    ymax = nodes(elements(e, 4),2);

    % Integrate stiffness using Gauss quadrature
    for g = 1 : size(xg, 1)
        % Map integration point from standard space (r, s) to cartesian coordinates (x, y)
        xgi(1) = 0.5 * (xmax - xmin) * xg(g, 1) + 0.5 * (xmax + xmin);
        xgi(2) = 0.5 * (ymax - ymin) * xg(g, 2) + 0.5 * (ymax + ymin);

        % Get B matrix for stiffness integration
        B = bMatrix(xmin, xmax, ymin, ymax, xgi);

        % Evaluate summand for quadrature of element stiffness matrix
        KE = KE + B' * C * B * wg(g);
    end
    KE = (xmax - xmin)/2 * (ymax - ymin)/2 * KE;

    % Use element connectivity to add element stiffness to global stiffness matrix
    indices = [2 * elements(e, :) - 1; 2 * elements(e, :)];
    indices = indices(:);
    K(indices, indices) = K(indices, indices) + KE;
end

% Load
f([10, 16]) = 1 / 2 * (700 * ly);

% Solve system
fixeddofs = [1, 2, 4];
freedofs = setdiff(1 : nNodes * dofPerNode, fixeddofs);

u(freedofs) = K(freedofs, freedofs) \ f(freedofs);

% Plot displacements
scale = 5;
uNodes = transpose(reshape(u, dofPerNode, nNodes));
plotElements(nodes + scale * uNodes, elements, 'Deformed system');

%% Postprocessing: Strain and stress
% Prepare plots
patchesSigma1 = plotElements(nodes, elements, 'Stress \sigma_1');
patchesSigma2 = plotElements(nodes, elements, 'Stress \sigma_2');
patchesSigma12 = plotElements(nodes, elements, 'Stress \tau');

% Initialize cells or arrays for stress
strainGauss = cell(nElements, 1);
stressGauss = cell(nElements, 1);
strainNodes = cell(nElements, 1);
stressNodes = cell(nElements, 1);
misesGauss = zeros(nElements, 4);
misesNodes = zeros(nElements, 4);

% Element-wise stress calculation
for e = 1 : nElements
    xmin = nodes(elements(e, 1), 1);
    xmax = nodes(elements(e, 3), 1);
    ymin = nodes(elements(e, 1), 2);
    ymax = nodes(elements(e, 3), 2);

    indices = [2 * elements(e, :) - 1; 2 * elements(e, :)];
    indices = indices(:);

    % Compute stress at integration points
    strainGauss{e} = zeros(4, 3);
    stressGauss{e} = zeros(4, 3);
    for g = 1 : size(xg, 1)
        % Map integration point from standard space (r, s) to cartesian coordinates (x, y)
        xgi(1) = 0.5 * (xmax - xmin) * xg(g, 1) + 0.5 * (xmax + xmin);
        xgi(2) = 0.5 * (ymax - ymin) * xg(g, 2) + 0.5 * (ymax + ymin);

        % Get B matrix for stress calculation
        B = bMatrix(xmin, xmax, ymin, ymax, xgi);

        % Evaluate stress at integration point
        strainGauss{e}(g, :) = B * u(indices);
        stressGauss{e}(g, :) = C * B * u(indices);

        % Compute von Mises stress at integration point
        misesGauss(e, g) = vonMisesStress(stressGauss{e}(g, :));
    end

    % Compute stress at nodes
    strainNodes{e} = zeros(4, 3);
    stressNodes{e} = zeros(4, 3);
    elementNodes = nodes(elements(e, :), :);
    for g = 1 : size(elementNodes, 1)
        % Get B matrix for stress calculation
        B = bMatrix(xmin, xmax, ymin, ymax, elementNodes(g, :));

        % Evaluate stress at nodes
        strainNodes{e}(g, :) = B * u(indices);
        stressNodes{e}(g, :) = C * B * u(indices);
    end
    % Compute von Mises stress at nodes in one go
    misesNodes(e, :) = vonMisesStress(stressNodes{e});

    patchesSigma1(e).set('FaceColor', 'interp', 'FaceVertexCData', stressNodes{e}(:, 1));
    patchesSigma1(e).FaceAlpha = 1;
    patchesSigma2(e).set('FaceColor', 'interp', 'FaceVertexCData', stressNodes{e}(:, 2));
    patchesSigma2(e).FaceAlpha = 1;
    patchesSigma12(e).set('FaceColor', 'interp', 'FaceVertexCData', stressNodes{e}(:, 3));
    patchesSigma12(e).FaceAlpha = 1;
end

% Averaging nodal stress
stressAvg = zeros(nNodes, 3);
stressAvg([1, 2], :) = stressNodes{1}([1 2], :);
stressAvg(3, :) = 1/2 * (stressNodes{1}(4, :) + stressNodes{2}(1, :));
stressAvg(4, :) = 1/3 * (stressNodes{1}(3, :) + stressNodes{2}(2, :) + stressNodes{3}(1, :));
stressAvg(5, :) = stressNodes{3}(2, :);
stressAvg(6, :) = stressNodes{2}(4, :);
stressAvg(7, :) = 0.5 * (stressNodes{2}(3, :) + stressNodes{3}(4, :));
stressAvg(8, :) = stressNodes{3}(3, :);

% Averaging von Mises stress
%   If we want to compute the von Mises stress at the nodes, the question about the order of computations
%   arises. We could compute the nodal stress for each element, average the nodal stresses and then compute the
%   von Mises stress from the averaged nodal stress. On the other hand, we could first compute the nodal von
%   Mises stress for each element and then average the von Mises stress, similar as for the nodes. This process
%   is shown in the following, as we computed the von Mises stress at the nodes already.
% Abaqus: Result options -> Averaging -> Compute scalars BEFORE averaging (default)
misesAvg = zeros(nNodes, 1);
misesAvg([1 2]) = misesNodes(1, [1 2]);
misesAvg(3) = 0.5 * (misesNodes(1, 4) + misesNodes(2, 1));
misesAvg(4) = 1/3 * (misesNodes(1, 3) + misesNodes(2, 2) + misesNodes(3, 1));
misesAvg(5) = misesNodes(3, 2);
misesAvg(6) = misesNodes(2, 4);
misesAvg(7) = 0.5 * (misesNodes(2, 3) + misesNodes(3, 4));
misesAvg(8) = misesNodes(3, 3);

%% Plots
plotElements(nodes, elements, 'Averaged \sigma_1', stressAvg(:, 1));

%   Two other functions for plotting are included:
%   plotField() interpolates the given values at nodes using the bilinear shape functions. The primary example
%       is the solution u(x) = sum of (h_I * u_I), but we can interpolate between any given nodal quantity in
%       the same way. The field function is used by Matlab's fsurf() function for plotting.
%
%       Thereby, feeding the averaged sigma_1 to plotElements() and to plotField() should produce similar
%       results, where in plotElements() the Matlab patch interpolates the color between nodal values while
%       plotField() uses our shape functions to explicitly interpolate the stress between the nodes.
plotField(nodes, elements, stressAvg(:, 1), 'Averaged \sigma_1 field');

%   plotStressField() evaluates the stress field at a given point for plotting using Matlab's fsurf().
plotStressField(nodes, elements, C, u, 1, '\sigma_1 field');

%% Additional functions
% Shape functions, H and B matrices
function h = shapeFunctions(xmin, xmax, ymin, ymax, x)
    % Evaluate shape functions for rectangular element at x in cartesian coordinates
    % INPUT
    %   xmin (scalar)       Minimum x coordinate for element
    %   xmax (scalar)       Maximum x coordinate for element
    %   ymin (scalar)       Minimum y coordinate for element
    %   ymax (scalar)       Maximum y coordinate for element
    %   x (array)           Point(s) for evaluation of shape functions in cartesian coordinates

    A = (xmax - xmin) * (ymax - ymin);

    % Evaluate shape functions
    h(:, 1) = 1/A * (xmax - x(:, 1)) .* (ymax - x(:, 2));
    h(:, 2) = 1/A * (x(:, 1) - xmin) .* (ymax - x(:, 2));
    h(:, 3) = 1/A * (x(:, 1) - xmin) .* (x(:, 2) - ymin);
    h(:, 4) = 1/A * (xmax - x(:, 1)) .* (x(:, 2) - ymin);
end

function dh = shapeFunctionGradients(xmin, xmax, ymin, ymax, x)
    % Evaluate shape functions gradients for rectangular element at x in cartesian coordinates
    % INPUT
    %   xmin (scalar)       Minimum x coordinate for element
    %   xmax (scalar)       Maximum x coordinate for element
    %   ymin (scalar)       Minimum y coordinate for element
    %   ymax (scalar)       Maximum y coordinate for element
    %   x (array)           Point(s) for evaluation of shape functions in cartesian coordinates

    A = (xmax - xmin) * (ymax - ymin);

    dh(1:2, 1) = -1/A * [(ymax - x(2)) ; (xmax - x(1))];
    dh(1:2, 2) = 1/A * [(ymax - x(2)); (-(x(1) - xmin))];
    dh(1:2, 3) = 1/A * [(x(2) - ymin); ((x(1) - xmin))];
    dh(1:2, 4) = -1/A * [(x(2) - ymin); -(xmax - x(1))];
end

function H = hMatrix(xmin, xmax, ymin, ymax, x)
    % Evaluate H matrix for linear quadrilaterial element
    % INPUT
    %   xmin (scalar)   Minimum x coordinate for element
    %   xmax (scalar)   Maximum x coordinate for element
    %   ymin (scalar)   Minimum y coordinate for element
    %   ymax (scalar)   Maximum y coordinate for element
    %   x (array)           Point(s) for evaluation of shape functions in cartesian coordinates
    % OUTPUT
    %   H (array)       H matrix for linear quadrilateral element

    H = [];
end

function B = bMatrix(xmin, xmax, ymin, ymax, x)
    % Evaluate H matrix for linear quadrilaterial element
    % INPUT
    %   xmin (scalar)   Minimum x coordinate for element
    %   xmax (scalar)   Maximum x coordinate for element
    %   ymin (scalar)   Minimum y coordinate for element
    %   ymax (scalar)   Maximum y coordinate for element
    %   x (array)           Point(s) for evaluation of shape functions in cartesian coordinates
    % OUTPUT
    %   B (array)       B matrix for linear quadrilateral element

    dh = shapeFunctionGradients(xmin, xmax, ymin, ymax, x);
    B = [dh(1,1), 0, dh(1,2), 0, dh(1,3), 0, dh(1,4), 0;
        0, dh(2,1), 0, dh(2,2), 0, dh(2,3), 0, dh(2,4);
        dh(2,1), dh(1,1), dh(2,2), dh(1,2), dh(2,3), dh(1,3), dh(2,4), dh(1,4)];
end

% Mises
function mises = vonMisesStress(stress)
    % Compute mises stress from given stress vector in voigt notation
    % INPUT
    %   stress (array)          Stress vector(s) as row vector
    % OUTPUT
    %   mises (scalar/array)    von Mises stress for each row of stress input array
    
    mises = sqrt(stress(1).^2 + stress(2).^2 - stress(1) * stress(2) + 3 * stress(3)^2);
end

% Field functions

% These are functions to evaluate quantities either by interpolating from nodal values or 
% directly calculating the stress component at a given position.

function u = field(xmin, xmax, ymin, ymax, x, y, nodalValues)
    % Interpolate nodal values to given location
    % INPUT
    %   xmin (scalar)   Minimum x coordinate for element
    %   xmax (scalar)   Maximum x coordinate for element
    %   ymin (scalar)   Minimum y coordinate for element
    %   ymax (scalar)   Maximum y coordinate for element
    %   x (array)       x coordinate(s) for evaluation of shape functions in cartesian coordinates
    %   y (array)       y coordinate(s) for evaluation of shape functions in cartesian coordinates
    % OUPUT
    %   u (scalar)      Interpolated quantity

    u = zeros(size(x));
    H = shapeFunctions(xmin, xmax, ymin, ymax, [x(:), y(:)]);
    for i = 1 : numel(u)
        u(i) = H(i, :) * nodalValues;
    end
    % u = reshape(u, size(x));
end

function s = stressField(xmin, xmax, ymin, ymax, x, y, C, u, component)
    % Compute stress tensor at location x and return specific component
    % INPUT
    %   xmin (scalar)       Minimum x coordinate for element
    %   xmax (scalar)       Maximum x coordinate for element
    %   ymin (scalar)       Minimum y coordinate for element
    %   ymax (scalar)       Maximum y coordinate for element
    %   x (array)           x coordinate(s) for evaluation of stress
    %   y (array)           y coordinate(s) for evaluation of stress
    %   C (array)           Stress-strain or constitutive matrix
    %   u (array)           Vector of nodal displacements
    %   component (scalar)  0: Mises, 1: sigma_x, 2: sigma_y, 3: sigma_xy
    % OUPUT
    %   s (scalar/array)    Stress value(s) for given location

    s = zeros(size(x));
    for i = 1 : numel(s)
        B = bMatrix(xmin, xmax, ymin, ymax, [x(i), y(i)]);
        stress = C * B * u;
        if component == 0
            % Mises
            s(i) = vonMisesStress(stress');
        else
            s(i) = stress(component);
        end
    end
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

function plotStressField(nodes, elements, C, u, component, title)
    % Plot stress field for given component
    % INPUT
    %   nodes (array)       Nodal coordinates as (nNodes x 2) array
    %   elements (array)    Nodal indices per element as (nElements x 4) array
    %   C (array)           Stress-strain or constitutive material matrix
    %   u (array)           Vector of nodal displacements of size (DOF x 1)
    %   component (scalar)  Stress component for plotting
    %   title (char)        Title for plot

    newFigure(title);
    hold on;
    for e = 1 : 3
        xmin = nodes(elements(e, 1), 1);
        xmax = nodes(elements(e, 3), 1);
        ymin = nodes(elements(e, 1), 2);
        ymax = nodes(elements(e, 3), 2);
        indices = [2 * elements(e, :) - 1; 2 * elements(e, :)];
        indices = indices(:);
        fhStress = @(x, y) stressField(xmin, xmax, ymin, ymax, x, y, C, u(indices), component);
        fsurf(fhStress, [xmin xmax ymin ymax], 'LineStyle', 'none');
    end
    hold off;
    view(2)
    xlim([min(nodes(:,1)), max(nodes(:,1))]);
    ylim([min(nodes(:,2)), max(nodes(:,2))]);
end

function plotField(nodes, elements, nodalValues, title)
    % Plot field by interpolating between nodal values
    % INPUT
    %   nodes (array)       Nodal coordinates as (nNodes x 2) array
    %   elements (array)    Nodal indices per element as (nElements x 4) array
    %   nodalValues (array) Nodal values for interpolation and plotting
    %   title (char)        Title for plot

    newFigure(title);
    hold on;
    for e = 1 : 3
        xmin = nodes(elements(e, 1), 1);
        xmax = nodes(elements(e, 3), 1);
        ymin = nodes(elements(e, 1), 2);
        ymax = nodes(elements(e, 3), 2);
        fhf = @(x, y) field(xmin, xmax, ymin, ymax, x, y, nodalValues(elements(e, :)));
        fsurf(fhf, [xmin xmax ymin ymax], 'LineStyle', 'none');
    end
    hold off;
    view(2)
    xlim([min(nodes(:,1)), max(nodes(:,1))]);
    ylim([min(nodes(:,2)), max(nodes(:,2))]);
    clim([min(nodalValues), max(nodalValues)])
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