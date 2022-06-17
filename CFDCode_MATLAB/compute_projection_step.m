
%% Appendix E: compute_projection_step.m

function [U_n1,V_n1] = compute_projection_step(U_star,V_star)
    
    global dx dy Nx Ny K1_prime K2_prime
    
    dudx = zeros(Nx,Ny);
    dudx(1:end-1,:) = (U_star(2:end,:)-U_star(1:end-1,:))./dx;
    dudx(end,:) = (U_star(1,:)-U_star(end,:))./dx;
    
    dvdy = zeros(Nx,Ny);
    dvdy(:,1:end-1) = (V_star(:,2:end)-V_star(:,1:end-1))./dy;
    dvdy(:,end) = (V_star(:,1)-V_star(:,end))./dy;
    
    f_hat = fftshift(fft2(dudx+dvdy)./(Nx.*Ny));
    P_hat = -f_hat./(K1_prime.^2+K2_prime.^2);
    P_hat(Nx./2+1,Ny./2+1) = 0;
    
    P = ifft2(fftshift(P_hat)).*(Nx.*Ny);
    
    dPdx = zeros(Nx,Ny);
    dPdx(2:end,:) = (P(2:end,:)-P(1:end-1,:))./dx;
    dPdx(1,:) = (P(1,:)-P(end,:))./dx;
    U_n1 = -dPdx + U_star;
    
    dPdy = zeros(Nx,Ny);
    dPdy(:,2:end) = (P(:,2:end)-P(:,1:end-1))./dy;
    dPdy(:,1) = (P(:,1)-P(:,end))./dy;
    V_n1 = -dPdy + V_star;
    
end
