function Utap=tapper(U, wx, wz, x0, z0, x, z,npml)
%% wx, wz is the width of tapper; x0, z0 is the start of tapper, usually at the locations of srcs and recs

    Utap = U;
    dx = x(2) - x(1);
    dz = z(2) - z(1);
    nx = length(x);
    nz = length(z);
    nwx = fix(wx / dx);
    nwz = fix(wz / dz);
    lx = 0:nwx-1;
    lz = 0:nwz-1;
    xcos = (1 - cos(lx/nwx * pi))/2;
    zcos = (1 - cos(lx/nwz * pi))/2;
    xcos = [zeros(1, npml+fix(x0/dx)), xcos];
    zcos = [zeros(1, npml+fix(z0/dz)), zcos];
    if wx ~= 0
    Utap(1:nwx + npml + fix(x0/dx),:) = Utap(1:nwx + npml + fix(x0/dx),:) .* repmat(xcos', 1, nz);
    Utap(end:-1:end-nwx+1-npml -  fix(x0/dx),:) = Utap(end:-1:end-nwx+1-npml -  fix(x0/dx),:) .* repmat(xcos', 1, nz);   
    end
    if wz ~= 0
    Utap(:,1:nwz + npml  + fix(z0/dz)) = Utap(:,1:nwz + npml  + fix(z0/dz)) .* repmat(xcos, nx, 1);
    Utap(:,end:-1:end-nwz+1 - npml -  fix(z0/dz)) = Utap(:,end:-1:end-nwz+1-npml- fix(z0/dz)) .* repmat(xcos, nx, 1);  
    end
end