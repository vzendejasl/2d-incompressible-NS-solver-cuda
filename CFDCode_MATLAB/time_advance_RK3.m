
%% Appendix A: time_advance_RK3.m

function [U,V] = time_advance_RK3(U0,V0,dt)
    
    [U_k1,V_k1] = compute_RHS(U0,V0);
    U_k1 = U_k1.*dt;
    V_k1 = V_k1.*dt;

    [U_k2,V_k2] = compute_projection_step(U0+U_k1./2,V0+V_k1./2);
    [U_k2,V_k2] = compute_RHS(U_k2,V_k2);
    U_k2 = U_k2.*dt;
    V_k2 = V_k2.*dt;
    
    [U_k3,V_k3] = compute_projection_step(U0-U_k1+2.*U_k2,V0-V_k1+2.*V_k1);
    [U_k3,V_k3] = compute_RHS(U_k3,V_k3);
    U_k3 = U_k3.*dt;
    V_k3 = V_k3.*dt;

    [U,V] = compute_projection_step(U0 + (U_k1 + 4.*U_k2 + U_k3)./6, V0 + (V_k1 + 4.*V_k2 + V_k3)./6);
    
end
