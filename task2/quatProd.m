function qprod = quatProd(ql, qr)

    if numel(ql) == 3 % assume pure quat
        ql = [0; ql];
    end
    
    if numel(qr) == 3 % assume pure quat
        qr = [0; qr];
    end

    etal  = ql(1); 
    epsl  = ql(2:4);
    etar  = qr(1); 
    epsr  = qr(2:4);

    qprod = [etal*etar-epsl'*epsr
             etar*epsl+etal*epsr+cross(epsl,epsr)];
end