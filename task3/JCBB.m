function abest = JCBB(z, zbar, S, alpha1, alpha2)
    m = size(z, 1)/2;
    abest = zeros(m, 1);
    a = zeros(m, 1);
    ic = individualCompatibility(z, zbar, S);
    g2 = chi2inv(1 - alpha2, 2);
    % reorder to remove problems with first meas being outside conf
    [~, order] = sort(min(ic, [], 2));
    j = 1;
    zo = z( reshape([2*order - 1, 2*order]',[], 1));
    ico = ic(order, :);
    
    abesto = JCBBrec(zo, zbar, S, alpha1, g2, j, a, ico, abest);
    abest(order) = abesto;
end

function abest = JCBBrec(z, zbar, S, alpha1, g2, j, a, ic, abest)
    m = size(z, 1)/2;
    n = nnz(a);
    
    if j > m % end of recursion, no more measurements to associate
        if n > nnz(abest) || ((n >= nnz(abest)) && (NIS(z, zbar, S, a) < NIS(z, zbar, S, abest)))
            abest = a;
        % else abest = previous abest from the input
        end
    else  % still at least one measurement to associate
        [~, I] = sort(ic(j, ic(j, :) < g2)); % be gready, indexing wrong...
        allinds = 1:size(ic, 2);
        usableinds = allinds(ic(j, :) < g2);
        for i = usableinds(I)
            a(j) = i;
            if NIS(z, zbar, S, a) < chi2inv(1-alpha1, 2 * (n+1)) % jointly compatible?
                ici = ic(j:end, i);
                ic(j:end, i) = inf; % landmark not available any more.
                abest = JCBBrec(z, zbar, S, alpha1, g2, j + 1, a, ic, abest);
                ic(j:end, i) = ici; % set landmark available again for next round.
            end
        end
        if n + (m - j - 1) >= nnz(abest) % try this measurement as clutter, should one use poisson assumptions here...?
            a(j) = 0;
            abest = JCBBrec(z, zbar, S, alpha1, g2, j + 1, a, ic, abest);
        end
    end
end

function ic = individualCompatibility(z, zbar, S)
    ic = zeros(size(z,1)/2, size(zbar,1)/2);
    for i = 1:size(zbar, 1)/2
        for j = 1:size(z, 1)/2
            v = z(2*(j-1) + [1;2])- zbar(2*(i-1) + [1;2]);
            %v(2) = wrapToPi(v(2));
            ic(j, i) = (v' / S(2*(i-1) + [1,2], 2*(i-1) + [1,2])) * v;
        end
    end
end

function nis = NIS(z, zbar, S, a)
    zr = reshape(z, 2, []);
    zbarr = reshape(zbar, 2, []);
    
    ztest = zr(:, a > 0);
    zbartest = zbarr(:, a(a> 0));
    
    inds = 2*a(a>0) - 1;
    inds = [inds, inds + 1]'; inds = inds(:);
    Stest = S(inds, inds);
    
    v = ztest - zbartest;
    v = v(:);
    v(2:2:end) = wrapToPi(v(2:2:end));
    
    nis = (v' / Stest) * v;
end