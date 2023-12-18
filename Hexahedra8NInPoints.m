function [weig, posgp, shapef, dershapef] = Hexahedra8NInPoints
    % Hexahedra8NInPoints: Shape function routine for Hexahedra8 element

    weig = ones(1, 8); % Using ones to create a vector of ones
    posgp = 1/sqrt(3) * [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1]';


    shapef = zeros(8, 8);
    dershapef = zeros(3, 8, 8);

    for j = 1:8
        xi = posgp(:, j)';
        [N, B] = shape_functions(xi);
        shapef(j, :) = N;
        dershapef(:, :, j) = B;
    end
end

function [N, B] = shape_functions(xii)
    % Hexahedra8N: Shape functions and derivatives for Hexahedra8 element

    xi = xii(1);
    eta = xii(2);
    zeta = xii(3);

    N = 0.125 * [(1-xi)*(1-eta)*(1-zeta), (1+xi)*(1-eta)*(1-zeta), (1+xi)*(1+eta)*(1-zeta), (1-xi)*(1+eta)*(1-zeta), (1-xi)*(1-eta)*(1+zeta), (1+xi)*(1-eta)*(1+zeta), (1+xi)*(1+eta)*(1+zeta), (1-xi)*(1+eta)*(1+zeta)];

    B = 0.125 * [-(1-eta)*(1-zeta), (1-eta)*(1-zeta), (1+eta)*(1-zeta), -(1+eta)*(1-zeta), -(1-eta)*(1+zeta), (1-eta)*(1+zeta), (1+eta)*(1+zeta), -(1+eta)*(1+zeta);
                -(1-xi)*(1-zeta), -(1+xi)*(1-zeta), (1+xi)*(1-zeta), (1-xi)*(1-zeta), -(1-xi)*(1+zeta), -(1+xi)*(1+zeta), (1+xi)*(1+zeta), (1-xi)*(1+zeta);
                -(1-xi)*(1-eta), -(1+xi)*(1-eta), -(1+xi)*(1+eta), -(1-xi)*(1+eta), (1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
end




