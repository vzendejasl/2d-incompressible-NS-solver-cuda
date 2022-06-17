
%% Appendix B: Order of Accuracy Test of RK3 Scheme

clear
close all

dt_vec = [1 5e-1 1e-1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4];
t_final = 1;
U0 = 5;
error = zeros(1,length(dt_vec));

for i = 1:length(dt_vec)
    
    dt = dt_vec(i);
    
    U_RK3 = U0;
    diff = zeros(1,t_final./dt+1);
    for t = dt:dt:t_final
        U_exact = U0.*exp(t);
        U_RK3 = time_advance_RK3_test(U_RK3,dt);
        diff = abs(U_exact-U_RK3);
    end
    
    error(i) = max(diff);
    
end

figure
plot(dt_vec,error,'-b','LineWidth',2);
grid on
xlabel('dt');
ylabel('max_{0<=t<=1}(|u_{exact} - u_{RK3}|)');
set(gca,'XScale','log','YScale','log');

RK3_order = (log(error(end))-log(error(1)))./(log(dt_vec(end))-log(dt_vec(1)));
disp("Slope of RK3 error = " + num2str(RK3_order));

function [U] = time_advance_RK3_test(U0,dt)
    
    U_k1 = U0.*dt;
    U_k2 = (U0+U_k1./2).*dt;
    U_k3 = (U0-U_k1+2.*U_k2).*dt;
    U = U0 + (U_k1 + 4.*U_k2 + U_k3)./6;

end
