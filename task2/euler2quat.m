function q = euler2quat(eul)% (phi, theta, psi)
    % finds the quaternion corresponding to the rotation roll pitch yaw
    % (phi, theta, psi) in radians.
    
    halfAngles = 0.5*eul;
    CS = [cos(halfAngles), sin(halfAngles)];

    hasSinFirst = diag(ones(3,1), -1);
    hasSinFirst = hasSinFirst(:, 1:(end-1)); % binary matrix if the first part of (10.39) has sin for an angle
    hasSinLast = (ones(4,3) - hasSinFirst); % last part has sin
    
    F = (1:3) + 3 * hasSinFirst; % indices into CS for needed numbers of the first part
    L = (1:3) + 3 * hasSinLast; % indices for second part
    sign = diag([1, -1, 1, -1]); % sign matrix for second part
    q = prod(CS(F), 2) + sign * prod(CS(L), 2); % perform product, get correct sign, and sum.
end