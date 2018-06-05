clc, clear
%% ================= Simulation of LODCO Algorithm =================

%% 基本参数设置
k = 1e-28;                        % 有效开关电容
tau = 0.002;                      % 时间片长度(s)
tau_d = 0.002;                    % 计算任务执行时间的deadline(s)
phi = 0.002;                      % 任务丢弃的惩罚项权重(s)
omega = 1e6;                      % 服务器带宽(Hz)
sigma = 1e-13;                    % 接收端的噪声功率(W)
p_tx_max = 1;                     % 移动设备最大传输功率(W)
f_max = 1.5e9;                    % 移动设备最大CPU时钟周期频率(Hz)
E_max = 0.002;                    % 电池允许的最大放电量(J)
L = 1000;                         % 一项计算任务的大小(bit)
X = 737.5;                        % 移动设备执行1bit任务所需的时钟周期个数
W = L*X;                          % 移动设备本地执行一项计算任务所需的时钟周期个数
E_H_max = 48e-6;                  % 收集的能量服从的均匀分布上限(J)
g0 = power(10,-4);                % 路径损失常数(dB转化之后的数值比)
d = 50;                           % 服务器和移动设备之间的相对距离(m)

%% 变量控制
T = 3000;                         % 时间片个数
E_min = 0.02e-3;                  % 能量使用下界(J)
V = 1e-5;                         % LODCO中penalty项的权重(J^2/second)
rho = 0.6;                        % 计算任务抵达的概率

% 实际能耗上限
E_max_hat = min(max(k*W*(f_max^2), p_tx_max*tau), E_max);
theta = E_max_hat + V*phi/E_min;        % 扰动参数


%% 中间变量存储
B = zeros(T, 1);                        % 实际电量
B_hat = zeros(T, 1);                    % 虚拟电量
e = zeros(T, 1);                        % 能量收集
indicator = zeros(T, 1);                % 1/2/3/4分别对应local,remote,dropped,无任务产生
f = zeros(T, 1);                        % 移动设备本地执行的频率
p = zeros(T, 1);                        % 移动设备卸载执行的传输功率
cost = zeros(T, 3);                     % 每一列分别对应local,remote,最后选择的模式下的execution cost
E = zeros(T, 3);                        % 每一列分别对应local,remote,最后选择的模式下的能耗

% 关闭fslove的输出
opt = optimset('Display', 'off');

t = 1;
while t <= T
    %% 阶段初始化
    % 以伯努利分布产生计算任务
    zeta = binornd(1, rho);
    % 产生虚拟电量队列值
    B_hat(t) = B(t) - theta;

    %% 求解optimal energy harvesting e* (不管是否有任务产生，能量收集都不能停！)
    % 产生E_H_t
    E_H_t = unifrnd(0, E_H_max);
    if B_hat(t) <= 0
        e(t) = E_H_t;
    end
    
    if zeta == 0
        % 没有计算任务产生
        indicator(t) = 4;
        % f(t) = 0; p(t) = 0;  默认即为0
    else
        % 产生信道功率增益
        h = exprnd(g0/power(d,4));

        %% 求解P_ME
        f_L = max(sqrt(E_min/(k*W)), W/tau_d);
        f_U = min(sqrt(E_max/(k*W)), f_max);
        if f_L <= f_U
            % P_ME有解
            f0 = power(V/(-2*B_hat(t)*k), 1/3);
            if f0 > f_U
                f(t) = f_U;
            elseif f0 >= f_L && f0 <= f_U && B_hat(t) < 0
                f(t) = f0;
            elseif f0 < f_L
                f(t) = f_L;
            end
            % 计算local的execution delay
            cost(t, 1) = W / f(t);
            % 计算此时的能耗
            E(t, 1) = k * W * (f(t)^2);
            if E(t, 1) >= B(t)
                disp(['本地执行电量不足![此时t为', num2str(t), ']']);
                % 当电量不足或者本地执行子问题无解的时候，将该子问题的目标函数值设为inf
                J_m = inf;
            else
                % 计算此时的J_m，即本地执行子问题的目标函数值
                J_m = -B_hat(t)*E(t,1) + V*cost(t,1);
            end

        else
            disp(['P_ME无解![此时t为', num2str(t), ']']);
            J_m = inf;
        end

        %% 求解P_SE
        E_tmp = sigma*L*log(2) / (omega*h);
        p_L_taud = (power(2, L/(omega*tau_d)) - 1) * (sigma/h);
        if E_tmp >= E_min
            p_L = p_L_taud;
        else
            % 计算p_Emin
            y = @(x)x*L-omega*log2(1+h*x/sigma)*E_min;
            %p_Emin = double(vpa(solve(y, 1)));
            tmp = fsolve(y, [0.001, 1], opt);
            p_Emin = max(tmp);
            p_L = max(p_L_taud, p_Emin);
        end
        if E_tmp >= E_max
            p_U = 0;
        else
            % 计算p_Emax
            y = @(x)x*L-omega*log2(1+h*x/sigma)*E_max;
            p_Emax = max(fsolve(y, [0.001, 100]));
            p_U = min(p_tx_max, p_Emax);
        end
        if p_L <= p_U
            % P_SE有解
            % 计算p0
            tmp = B_hat(t);
            syms x
            y = tmp*log2(1+h*x/sigma) + h*(V-tmp*x)/(log(2)*(sigma+h*x));
            p0 = double(vpasolve(y));
            if p_U < p0
                p(t) = p_U;
            elseif p_L > p0 && B_hat(t) < 0
                p(t) = p_L;
            elseif p_L <= p0 && p_U >= p0 && B_hat(t) < 0
                p(t) = p0;
            end
             % 计算achievable rate
            r = calAchieveRate(h, p(t), omega, sigma);
            % 计算此时的execution delay
            cost(t, 2) = L / r;
            % 计算此时的能耗
            E(t, 2) = p(t) * L / r;
            if E(t, 2) >= B(t)
                disp(['卸载执行电量不足![此时t为', num2str(t), ']']);
                J_s = inf;
            else
                % 计算此时的J_s卸载执行子问题的目标函数值
                J_s = -B_hat(t)*E(t,2) + V*cost(t,2);
            end
        else
            disp(['P_SE无解![此时t为', num2str(t), ']']);
            J_s = inf;
        end

        %% 选取最佳模式
        [~, mode] = min([J_m, J_s, phi]);
        indicator(t) = mode;
    end

    % 计算实际选择的模式下的execution cost和能耗
    if indicator(t) == 1
        cost(t, 3) = cost(t, 1);
        E(t, 3) = E(t, 1);
    elseif indicator(t) == 2
        cost(t, 3) = cost(t, 2);
        E(t, 3) = E(t, 2);
    % 剩下两种情况能耗为0，不变
    elseif indicator(t) == 3
        cost(t, 3) = phi;
    else
        % 即没有任务产生的情况
        cost(t, 3) = 0;
    end
    % 电量迭代
    B(t+1) = B(t) - E(t, 3) + e(t);
    % 时间片迭代
    t = t + 1;
end

%% results
% 产生任务的时间片个数
num = T - sum(indicator == 4);
disp(['当前任务到达的频率为：', num2str(num/T)]);
mode1 = sum(indicator==1) / (num);
disp(['任务本地执行的比率：', num2str(mode1)]);
mode2 = sum(indicator==2) / (num);
disp(['任务卸载执行的比率:', num2str(mode2)]);
mode3 = sum(indicator==3) / (num);
disp(['任务被抛弃的比率：', num2str(mode3)]);