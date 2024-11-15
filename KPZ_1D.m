clear;
clc;
format long
tic;

L_all = 20:10:60;
lambda = 1;
dx = 1;

t_max = 1e3;
dt = 0.025;
t = 0:dt:t_max;

T = 1;

M_measure = 1e3;

order_sqrt_final = zeros(length(L_all),length(t));

for n = 1:length(L_all)
    L = L_all(n);
    x = 0:dx:L;
    h = zeros(length(x),length(t));
    order_sqrt_final(n,:) = time_evoluition(h, t, M_measure, dt, dx, T, lambda);
end
toc;

order_sqrt_final_mean = mean(order_sqrt_final(:,round(3/5*length(t):end)),2);

le = cell(1, length(L_all));
for i = 1:length(L_all)
    le{i} = strcat('L = ', num2str(L_all(i)));
end

figure;
plot(log(t), log(order_sqrt_final))
legend(le)

figure;
plot(log(L_all), log(order_sqrt_final_mean))

save(strcat('KPZ_data_CH3_num', num2str(M_measure),'_lambda', num2str(lambda), '.mat'), '-v7.3');


function y = Heun(x, dt, dx, T, lambda)
    noise = randn(length(x),1) * sqrt(2 * T / (dt * dx));
    yt = x + f(x, dx) * dt + lambda * g_CH3(x, dx) * dt + noise * dt;
    y = x + (f(x, dx) + f(yt, dx)) * dt / 2 + lambda * (g_CH3(x, dx) + g_CH3(yt, dx)) * dt / 2 + noise * dt;
end

function y = f(x, dx)
    y = (circshift(x, 1) + circshift(x, -1) - 2 * x) / (dx)^2;
end

function y = g_basic(x, dx)
    y = ((circshift(x, 1) - circshift(x, -1)) / (2 * dx)) .^ 2;
end

function y = g_CH3(x, dx)
    y = (((circshift(x, 1) - x) / dx) .^ 2 + ((circshift(x, -1) - x) / dx) .^ 2) / 2;
end

function y = g_LS(x, dx)
    x_l = circshift(x, 1) - x;
    x_r = x - circshift(x, -1);
    y = (x_l .^ 2 + x_r .^ 2 + x_l.*x_r) / (3*dx^2);
end


function order_final_sqrt = time_evoluition(h, t, M_measure, dt, dx, T, lambda)
    order_final = zeros(1,length(t));
    for n = 1:M_measure
        for i = 1:length(t)-1
            h(:,i+1) = Heun(h(:,i), dt, dx, T, lambda);
        end
        order_mean = mean(h);
        order2mean = mean(h.^2);
        order = order2mean - order_mean.^2;
        order_final = order_final + order;
    end
    order_final_sqrt = sqrt(order_final)/ M_measure;
end