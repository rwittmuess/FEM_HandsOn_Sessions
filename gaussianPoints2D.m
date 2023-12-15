function [xg, wg] = gaussianPoints2D(n)
    % Return integration points and weights for 2D gauss integration of order n
    % INPUT
    %   n (scalar)      OPTIONAL, order for integration, default: 2
    % OUTPUT
    %   xg (array)      Gaussian points in 2D
    %   wg (array)      Combined weights for 2D integration points

    if nargin < 1
        n = 2;
    end
    
    switch n
        case 1
            xg = [0, 0];
            wg = 2 * 2;
        case 2
            xg = [-1/sqrt(3), -1/sqrt(3);
                   1/sqrt(3), -1/sqrt(3);
                   1/sqrt(3),   1/sqrt(3);
                  -1/sqrt(3),  1/sqrt(3);];
            wg = ones(4,1);
        case 3
            xg = [-sqrt(3/5), -sqrt(3/5);
                  0, -sqrt(3/5);
                  sqrt(3/5), -sqrt(3/5);
                  -sqrt(3/5), 0;
                  0, 0;
                  sqrt(3/5), 0;
                  -sqrt(3/5), sqrt(3/5);
                  0, sqrt(3/5);
                  sqrt(3/5), sqrt(3/5)];
            wg = [5/9 * 5/9;
                  8/9 * 5/9;
                  5/9 * 5/9;
                  5/9 * 8/9;
                  8/9 * 8/9;
                  5/9 * 8/9;
                  5/9 * 5/9;
                  8/9 * 5/9;
                  5/9 * 5/9];
        otherwise
            error('Integration points for order %d not implemented', n);
    end
end