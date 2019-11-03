function R = quat2rotmat(quat)  
    thresh = 1e-6;
    if abs(norm(quat)-1)>thresh; error('input must be a unit quaternion. Got norm %f', norm(quat)); end

    eta = quat(1);
    eps = quat(2:4);

    S = crossProdMat(eps);
    R = eye(3) + 2*eta*S + 2*S^2;
end
