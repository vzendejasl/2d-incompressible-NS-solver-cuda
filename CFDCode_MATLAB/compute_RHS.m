
%% Appendix C: compute_RHS.m

function [U_RHS,V_RHS] = compute_RHS(U,V)

    global dx dy Re
    
    % U(i+1,j)
    U_ip1 = U;
    U_ip1(1:end-1,:) = U(2:end,:);
    U_ip1(end,:) = U(1,:);
    
    % U(i-1,j)
    U_im1 = U;
    U_im1(2:end,:) = U(1:end-1,:);
    U_im1(1,:) = U(end,:);
    
    % U(i,j+1)
    U_jp1 = U;
    U_jp1(:,1:end-1) = U(:,2:end);
    U_jp1(:,end) = U(:,1);
    
    % U(i,j-1)
    U_jm1 = U;
    U_jm1(:,2:end) = U(:,1:end-1);
    U_jm1(:,1) = U(:,end);
    
    % U(i+1,j-1)
    U_ip1_jm1 = U_jm1;
    U_ip1_jm1(1:end-1,:) = U_jm1(2:end,:);
    U_ip1_jm1(end,:) = U_jm1(1,:);
    
    % V(i+1,j)
    V_ip1 = V;
    V_ip1(1:end-1,:) = V(2:end,:);
    V_ip1(end,:) = V(1,:);
    
    % V(i-1,j)
    V_im1 = V;
    V_im1(2:end,:) = V(1:end-1,:);
    V_im1(1,:) = V(end,:);
    
    % V(i,j+1)
    V_jp1 = V;
    V_jp1(:,1:end-1) = V(:,2:end);
    V_jp1(:,end) = V(:,1);
    
    % V(i,j-1)
    V_jm1 = V;
    V_jm1(:,2:end) = V(:,1:end-1);
    V_jm1(:,1) = V(:,end);
    
    % V(i-1,j+1)
    V_im1_jp1 = V_im1;
    V_im1_jp1(:,1:end-1) = V_im1(:,2:end);
    V_im1_jp1(:,end) = V_im1(:,1);
    
    % U-Momentum
    duv_dy = ((V_jp1+V_im1_jp1).*(U_jp1+U)-(V+V_im1).*(U+U_jm1))./(4.*dy);
    duu_dx = (U_ip1.^2-U_im1.^2)./(2.*dx);
    d2u_dx2 = (U_ip1-2.*U+U_im1)./(dx.^2);
    d2u_dy2 = (U_jp1-2.*U+U_jm1)./(dy.^2);
    U_RHS = -duv_dy-duu_dx+(1./Re).*(d2u_dx2+d2u_dy2);
    
    % V-Momentum
    duv_dx = ((U_ip1+U_ip1_jm1).*(V_ip1+V)-(U+U_jm1).*(V+V_im1))./(4.*dx);
    dvv_dy = (V_jp1.^2-V_jm1.^2)./(2.*dy);
    d2v_dx2 = (V_ip1-2.*V+V_im1)./(dx.^2);
    d2v_dy2 = (V_jp1-2.*V+V_jm1)./(dy.^2);
    V_RHS = -duv_dx-dvv_dy+(1./Re).*(d2v_dx2+d2v_dy2);
    
end