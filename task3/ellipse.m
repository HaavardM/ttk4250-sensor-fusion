function ell = ellipse(mu, P, s, n)
thetas = (0:n)*2*pi/n;
ell = mu + s * chol(P)' * [cos(thetas); sin(thetas)];
end