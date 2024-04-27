function [u, v] = hs_estimation(Ix, Iy, It, tol, maxnum)
    % Ix,Iy,It: represent the spatial and temporal derivatives of intensity
    lambda_horn = 20; %regularization parameter
    [r, c] = size(Ix);

    % Defining the boundary condition matrix
    chorizontal = [3, ones(1, c-2)*8, 3];
    cvertical = [3; ones(r-2,1)*8; 3];
    cmtx = 8 * ones(r,c); 
    cmtx(1,:) = chorizontal;
    cmtx(r,:) = chorizontal;
    cmtx(:,1) = cvertical;
    cmtx(:,c) = cvertical;

    % Computing the elements needed for the update
    denominator = cmtx .* (Ix.^2 + Iy.^2) + lambda_horn * cmtx.^2;
    uv = (Ix .* Iy) ./ denominator;
    u1 = (Iy.^2 + lambda_horn * cmtx) ./ denominator;
    u2 = (Ix .* It) ./ denominator;
    v1 = (Ix.^2 + lambda_horn * cmtx) ./ denominator;
    v2 = (Iy .* It) ./ denominator;

    % Initializing the approximate optical flow
    u = zeros(r, c);
    v = zeros(r, c);

    k = 0;
    tot_err = 10^(7);

    % Main loop for iterative optimization
    while tot_err > tol && k < maxnum
        % Computing the temporary variables for the update
        %The coefficients are written for central difference scheme

        tmpu = conv2(u, [0 1 0; 1 1 1; 0 1 0], 'same');
        tmpv = conv2(v, [0 1 0; 1 1 1; 0 1 0], 'same');

        % Updating the flow estimates
        unew = u1 .* tmpu - uv .* tmpv - u2;
        vnew = v1 .* tmpv - uv .* tmpu - v2;

        % Computing the error and update the flow estimates
        tot_err = norm(unew - u, 'fro') + norm(vnew - v, 'fro');
        u = unew;
        v = vnew;
        k = k + 1;
    end
end
