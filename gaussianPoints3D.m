function [xg, wg] = gaussianPoints3D()
    % Return integration points and weights for 2D gauss integration of order 2
    % OUTPUT
    %   xg (array)      Gaussian points in 2D for full second order integration
    %   wg (array)      Combined weights for 2D integration points
    xg = [-1/sqrt(3), -1/sqrt(3), -1/sqrt(3);
           1/sqrt(3), -1/sqrt(3), -1/sqrt(3);
           1/sqrt(3),   1/sqrt(3), -1/sqrt(3);
          -1/sqrt(3),  1/sqrt(3), -1/sqrt(3);
          -1/sqrt(3), -1/sqrt(3), 1/sqrt(3);
           1/sqrt(3), -1/sqrt(3), 1/sqrt(3);
           1/sqrt(3),   1/sqrt(3), 1/sqrt(3);
          -1/sqrt(3),  1/sqrt(3), 1/sqrt(3)];
    wg = ones(8,1);
end