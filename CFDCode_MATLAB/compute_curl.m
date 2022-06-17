
%% Appendix G: compute_curl.m

function [curl] = compute_curl(U,V)
    
    global dx dy Nx Ny

    dvdx = zeros(Nx,Ny);
    dvdx(1:end-1,:) = (V(2:end,:)-V(1:end-1,:))./dx;
    dvdx(end,:) = (V(1,:)-V(end,:))./dx;
    
    dudy = zeros(Nx,Ny);
    dudy(:,1:end-1) = (U(:,2:end)-U(:,1:end-1))./dy;
    dudy(:,end) = (U(:,1)-U(:,end))./dy;
    
    curl = dvdx - dudy;
    
end

