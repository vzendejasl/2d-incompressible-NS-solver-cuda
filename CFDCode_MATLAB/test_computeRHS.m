
%% Appendix D: Order of Accuracy Test for 2D Spatial Solver

clear
close all

global dx dy Nx Ny Re

Nx_vec = [10 25 50 75 100 250 500 750 1000];
 
error_umom = zeros(1,length(Nx_vec));
error_vmom = zeros(1,length(Nx_vec));

for i = 1:length(Nx_vec)
    Nx = Nx_vec(i);
    
    Lx = 2*pi;
    Ny = Nx;
    Ly = Lx;
    Re = 1;

    dx = Lx/Nx;
    dy = Ly/Ny;

    % Define grid and exact solution for U-momentum equation
    x_umom = 0:dx:Lx-dx;
    y_umom = dy/2:dy:Ly-dy/2;
    [X_umom,Y_umom] = meshgrid(x_umom,y_umom);
    X_umom = X_umom';
    Y_umom = Y_umom';
    U_umom = cos(X_umom);
    U_RHS_exact = 2.*cos(X_umom).*sin(X_umom)-cos(X_umom).*cos(Y_umom)-1./Re.*cos(X_umom);

    % Define grid and exact solution for V-momentum equation
    x_vmom = dx/2:dx:Lx-dx/2;
    y_vmom = 0:dy:Ly-dy;
    [X_vmom,Y_vmom] = meshgrid(x_vmom,y_vmom);
    X_vmom = X_vmom';
    Y_vmom = Y_vmom';
    V_vmom = sin(Y_vmom);
    V_RHS_exact = -2.*cos(Y_vmom).*sin(Y_vmom)+sin(X_vmom).*sin(Y_vmom)-1./Re.*sin(Y_vmom);

    % Compute RHS for U-mom. and V-mom. equations using spatial solver
    [U_RHS_approx,V_RHS_approx] = compute_RHS(U_umom,V_vmom);
    
    error_umom(i) = max(max(abs(U_RHS_exact - U_RHS_approx)));
    error_vmom(i) = max(max(abs(V_RHS_exact - V_RHS_approx)));
    
end

figure
hold on
plot(Nx_vec,error_umom,'-b','LineWidth',2);
plot(Nx_vec,error_vmom,'--r','LineWidth',3);
hold off
grid on
xlabel('Number of grid points, N');
ylabel('max_{x,y}(|RHS_{exact} - RHS_{central FD}|)');
legend('U-momentum','V-momentum');
set(gca,'XScale','log','YScale','log');

Umom_order = (log(error_umom(end))-log(error_umom(1)))./(log(Nx_vec(end))-log(Nx_vec(1)));
disp("Slope of U-Momentum Error = " + num2str(Umom_order));
Vmom_order = (log(error_vmom(end))-log(error_vmom(1)))./(log(Nx_vec(end))-log(Nx_vec(1)));
disp("Slope of V-Momentum Error = " + num2str(Vmom_order));

