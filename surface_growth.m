clear;
clc;
format long
tic;

step = 2e4;
L = 100;
M_measure = 1;

order_mean = zeros(step+1,1);

for n = 1:M_measure

    S = zeros(L,step+1);
    order = zeros(step+1,1);
    diff = zeros(L-1,step+1);

    for i = 1:step
        S(:,i+1) = S(:,i);
        pos = ceil(rand*(L-2))+1;
        S(pos,i+1) = min(S(pos-1,i),S(pos+1,i))+1;
        % order(i+1) = sqrt(mean(S(:,i+1).^2)-mean(S(:,i+1)).^2);
        % order(i+1) = sqrt(mean(S(:,i+1).^2));
        order(i+1) = std(S(:,i+1));
        diff(:,i+1) = S(1:end-1,i+1) - S(2:end,i+1);
    end
    order_mean = order_mean + order;

end
order_mean = order_mean/M_measure;

figure;
plot(log(1:step+1),log(order_mean))

toc;