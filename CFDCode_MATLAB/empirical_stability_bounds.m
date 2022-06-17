
%% Appendix H: Empirical Stability Bounds

clear
close all

global dx dy Nx Ny Re K1_prime K2_prime

Nx = 64;
Lx = 2*pi;
Ny = Nx;
Ly = Lx;

dx = Lx/Nx;
dy = Ly/Ny;

% Staggered grid coordinates for U-momentum equation
x_umom = linspace(0,Lx-dx,Nx);
y_umom = linspace(dy/2,Ly-dy/2,Ny);
[X_umom,Y_umom] = meshgrid(x_umom,y_umom);
X_umom = X_umom';
Y_umom = Y_umom';

% Staggered grid coordinates for V-momentum equation
x_vmom = linspace(dx/2,Lx-dx/2,Nx);
y_vmom = linspace(0,Ly-dy,Ny);
[X_vmom,Y_vmom] = meshgrid(x_vmom,y_vmom);
X_vmom = X_vmom';
Y_vmom = Y_vmom';

% Modified wavenumbers
K1 = 2.*pi./Lx.*(-Nx/2:Nx/2-1).*ones(Nx);
K2 = (2.*pi./Ly.*(-Ny/2:Ny/2-1).*ones(Ny))';
K1_prime = 1./dx.*sqrt(2-2.*cos(K1.*dx));
K2_prime = 1./dy.*sqrt(2-2.*cos(K2.*dy));

% Initial Condition
U0 = rand(Nx,Ny);
V0 = rand(Nx,Ny);

Re_vec = [1 5 10 50 100 500 1000 5000 10000];
dt_vec = [1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 1e-1 5e-1 1];
stable = ones(length(Re_vec),length(dt_vec));

for Re_index = 1:length(Re_vec)
    
    Re = Re_vec(Re_index);
    
    for dt_index = 1:length(dt_vec)
        
        dt = dt_vec(dt_index);
        t = dt;
        t_final = dt*100;
        
        % Time loop
        U_num = U0;
        V_num = V0;
        while t <= t_final+dt
    
            % Numerical solution
            [U_num,V_num] = time_advance_RK3(U_num,V_num,dt);
            curl_num = compute_curl(U_num,V_num);
            
            if (max(max(abs(U_num-U0))) > 1e2) || (max(max(abs(V_num-V_num))) > 1e2)
                stable(Re_index,dt_index) = 0;
                break;
            end

            % Update time
            t = t + dt;
        end
    end
end
