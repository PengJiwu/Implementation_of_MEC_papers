clc, clear
%% ================= Simulation of LODCO-Based Greedy Algorithm =================

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
d0 = 1;                           % 服务器和移动设备之间的相对距离(m)

%% 变量控制
N = 20;                           % 移动设备数目
M = 4;                            % MEC服务器个数
T = 1000;                         % 时间片个数
E_min = 0.02e-3;                  % 能量使用下界(J)
V = 1e-4;                         % LODCO中penalty项的权重(J^2/second)
rho = 0.6;                        % 计算任务抵达的概率

% 实际能耗上限
E_max_hat = min(max(k*W*(f_max^2), p_tx_max*tau), E_max);
theta = E_max_hat + V*phi/E_min;        % 扰动参数

%% 中间变量存储           
B = zeros(T, N);                        % N个移动设备的实际电量
B_hat = zeros(T, N);                    % N个移动设备的虚拟电量
e = zeros(T, N);                        % N个移动设备的能量收集
indicator = zeros(T, N);                % 用1,2,3,4分别表示本地执行、卸载执行、drop(前三个均意味着有任务到达)以及没有任务到达
f = zeros(T, N);                        % N个移动设备本地执行的频率
p = zeros(T, N);                        % N个移动设备卸载执行的传输功率
local_execution_delay = zeros(T, N);    % N个移动设备的local execution delay
remote_execution_delay = zeros(T, N);   % N个移动设备的remote execution delay
cost = zeros(T, N);                     % N个移动设备最后选择的模式下的execution cost
E_local = zeros(T, N);                  % N个移动设备本地执行的能耗
E_remote = zeros(T, N);                 % N个移动设备卸载执行的能耗
E = zeros(T, N);                        % N个移动设备在最后选择的模式下的实际能耗
mode_num = zeros(T,3);                  % 每一列分别表示每轮中本地执行、远程执行及任务丢弃的比率(分母为有任务抵达的设备数)

% 关闭fslove的输出
opt = optimset('Display', 'off');

t = 1;
while t <= T
    % 用一个不定长的、3列的矩阵描述"移动设备i-MEC服务器j-i到j的最小卸载执行子问题目标函数值"
    map = [];
    % 存储每一个MEC服务器连接到的移动设备个数
    flags = zeros(M, 1);
    % 分别保存当前移动设备到各MEC服务器的J_s值
    J_s_matrix = zeros(N, M);
    % 分别存储每一个移动设备的本地执行的J_m，供二次决策时使用
    J_m = zeros(N, 1);
    % 产生虚拟电量队列值
    B_hat(t,:) = B(t,:) - theta;

    for i = 1:N
        %% 对每一个移动设备，阶段初始化
        %% 求解optimal energy harvesting e* (不管是否有任务产生，能量收集都不能停！)
        % 产生E_H_t
        E_H_t = unifrnd(0, E_H_max);
        if B_hat(t, i) <= 0
            e(t, i) = E_H_t;
        end

        % 以伯努利分布产生计算任务
        zeta = binornd(1, rho);
        if zeta == 0
            % 没有计算任务产生
            indicator(t, i) = 4;
            % f(t, i) = 0; p(t, i) = 0;    默认即为0
        else
            %% 求解P_ME
            f_L = max(sqrt(E_min/(k*W)), W/tau_d);
            f_U = min(sqrt(E_min/(k*W)), f_max);
            if f_L <= f_U
                % P_ME有解
                f0 = power(V/(-2*B_hat(t,i)*k), 1/3);
                if f0 > f_U
                    f(t, i) = f_U;
                elseif f0 >= f_L && f0 <= f_U && B_hat(t,i) < 0
                    f(t, i) = f0;
                elseif f0 < f_L
                    f(t, i) = f_L;
                end
                % 计算此时的execution delay
                local_execution_delay(t, i) = W / f(t, i);
                % 计算此时的能耗
                E_local(t, i) = k * W * (f(t, i)^2);
                if E_local(t, i) >= B(t, i)
                    disp(['本地执行电量不足![此时t为', num2str(t), ']']);
                    J_m(i) = inf;    % 设置为inf可以保证一定比phi小
                else
                    % 计算此时的J_m
                    % J_m代表的是子问题的目标函数值！是能耗和延迟统一变换后的结果，而非仅仅是延迟！
                    J_m(i) = -B_hat(t,i)*E_local(t, i) + V*local_execution_delay(t, i);
                end

            else
                disp(['P_ME无解![此时t为', num2str(t), ']']);
                J_m(i) = inf;
            end

            %% 求解P_SE
            % 随机产生服务器和移动设备的距离(限定在0 ~ 80之内)
            D = unifrnd(0, 100, N, M);
            % 服从lambda=1的指数分布的小尺度衰落信道功率收益
            gamma = exprnd(1, N, M);
            % 从任意移动设备到任意服务器的信道功率增益
            h = g0*gamma.*power(d0./D, 4);
            % 保存当前移动设备对所有服务器的卸载执行延迟
            remote_execution_delay_tmp = zeros(M, 1);

            for j = 1:M
                tmp_h = h(i,j);
                E_tmp = sigma*L*log(2) / (omega*tmp_h);
                p_L_taud = (power(2, L/(omega*tau_d)) - 1) * (sigma/tmp_h);
                if E_tmp >= E_min
                    p_L = p_L_taud;
                else
                    % 计算p_Emin
                    y = @(x)x*L-omega*log2(1+tmp_h*x/sigma)*E_min;
                    %p_Emin = double(vpa(solve(y, 1)));
                    tmp = fsolve(y, [0.0000000001, 1], opt);
                    p_Emin = real(max(tmp));
                    p_L = max(p_L_taud, p_Emin);
                end
                if E_tmp >= E_max
                    p_U = 0;
                else
                    % 计算p_Emax
                    %{
                    y = @(x)x*L-omega*log2(1+tmp_h*x/sigma)*E_max;
                    p_Emax = max(fsolve(y, [0.001, 100]));
                    p_U = min(p_tx_max, p_Emax);
                    %}
                    % 加速运算
                    p_Emax = 25;
                    p_U = 1;
                end
                if p_L <= p_U
                    % P_SE有解
                    % 计算p0
                    tmp = B_hat(t,i);
                    y = @(x)tmp*log2(1+tmp_h*x/sigma) + tmp_h*(V-tmp*x)/(log(2)*(sigma+tmp_h*x));
                    p0 = real(max(fsolve(y, [0.00000000001, 1], opt)));
                    if p_U < p0
                        p(t,i) = p_U;
                    elseif p_L > p0 && B_hat(t,i) < 0
                        p(t,i) = p_L;
                    elseif p_L <= p0 && p_U >= p0 && B_hat(t,i) < 0
                        p(t,i) = p0;
                    end
                    % 计算achievable rate
                    r = calAchieveRate(tmp_h, p(t,i), omega, sigma);
                    % 计算此时的延迟
                    remote_execution_delay_tmp(j) = L / r;
                    % 计算此时的能耗
                    E_remote(t, i) = p(t,i)*L / r;
                    if E_remote(t, i) >= B(t, i)
                        disp(['卸载执行电量不足![此时t为', num2str(t), ']']);
                        J_s = inf;
                    else
                        % 计算此时的J_s
                        % J_s是卸载执行子问题的目标函数值，而非延迟！！
                        J_s = -B_hat(t,i)*E_remote(t,i) + V*remote_execution_delay_tmp(j);
                    end
                else
                    disp(['P_SE无解![此时t为', num2str(t), ']']);
                    % 因此indicator(t, 2) = 0不变
                    J_s = inf;
                end
                J_s_matrix(i,j) = J_s;
            end
            
            % 选取出最佳的键值对组合
            [J_s_best, j_best] = min(J_s_matrix(i,:));
            % 计算此时的execution delay
            remote_execution_delay(t, i) = remote_execution_delay_tmp(j_best);

            %% 为第i个移动设备选取最佳模式
            [~, mode] = min([J_m(i), J_s_best, phi]);
            indicator(t, i) = mode;
            if mode == 2
                %{
                map(i) = j_best;
                map_Js_best = [map_Js_best, J_s_best];
                %}
                map = [map;[i,j_best,J_s_best]];
            end
        end
    end

    %% 为那些选择卸载执行的任务分配服务器(或另选模式)
    % UB = f_server_max*tau_d/(L*X),取f_server_max=f_max则约为4.0678
    UB = 4;
    while ~isempty(map)
        % 找到拥有最小传输延迟的移动设备-MEC服务器对
        [min_Js,index] = min(map(:,3));
        % 找到最小延迟对应的min_i和min_j
        min_i = map(index,1);
        min_j = map(index,2);

        % 判断此时对于移动设备min_i而言，卸载执行是否仍是最优mode
        [~, mode] = min([J_m(min_i), min_Js, phi]);
        % 当前最优模式仍为卸载执行
        if mode == 2
            % 当前移动设备很幸运，对其而言最优的卸载对象还可以再接收任务
            if flags(min_j) <= UB
                % 从map中删除该键值对并同步一系列共同维护的变量
                map(index,:) = [];
                % 对应的MEC服务器自增1(该操作的位置发生了变化，论文此处需要修改)
                flags(min_j) = flags(min_j) + 1;
            else
                % flags(min_j)已达上限，将那些最佳卸载对象为min_j的移动设备此时的目标函数值都设为inf，即不可达
                % (这个操作对于已经完成卸载的服务器而言不会有任何影响，因为它们不可能再出现在map中)
                % 此时这个设备min_i只能重新选择次最优的卸载对象
                J_s_matrix(:, min_j) = inf;
                % 此时无需更新map，之后map会根据最新的J_s_matrix来更新自己的值

                % 为min_i重新选择一个次最优的卸载对象.当该移动设备还有能选的服务器的时候，在这些可选的服务器中找到J_s最小的那个
                % 第一种情况: 还有服务器没满，可以选
                if min(J_s_matrix(min_i,:)) ~= inf
                    % 找到最小的那个之后，覆盖map中该移动设备选取的服务器，并覆盖对应的Js最小值
                    [min_Js_second, min_j_second] = min(J_s_matrix(min_i,:));
                    map(index,2:3) = [min_j_second,min_Js_second];
                    % 返回最外层的while，重新开始找当前具有最小J_s值的设备和对应的服务器
                    continue;
                % 第二种情况: 很不幸，每一个服务器都满了
                else
                    % 没有服务器可以选了，只能在另外两种mode中选取
                    % 重新设置指示变量
                    [~, mode] = min([J_m(min_i), inf, phi]);
                    indicator(t,i) = mode;
                    % 从map中删除该键值对并同步一系列共同维护的变量，因为map只存放欲图卸载执行的设备
                    map(index,:) = [];
                end
            end
        else
            % 设置新的最优模式，并从map中删除本键值对及维护的变量，因为map只存放欲图卸载执行的设备
            indicator(t,min_i) = mode;
            map(index,:) = [];
        end
    end

    % 计算每一个移动设备在当前模式下的execution cost
    cost(t,indicator(t,:)==1) = local_execution_delay(t,indicator(t,:)==1);
    cost(t,indicator(t,:)==2) = remote_execution_delay(t,indicator(t,:)==2);
    cost(t,indicator(t,:)==3) = phi;
    % 计算每一个移动设备在当前模式下的能耗
    E(t,indicator(t,:)==1) = E_local(t,indicator(t,:)==1);
    E(t,indicator(t,:)==2) = E_remote(t,indicator(t,:)==2);
    E(t,indicator(t,:)==3) = 0;

    % 中间结果输出
    task_num = N - size(find(indicator(t,:)==4),2);
    if task_num == 0
        local_rate = 0; offloading_rate = 0; drop_rate = 0;
    else
        local_rate = size(find(indicator(t,:)==1),2)/task_num;
        offloading_rate = size(find(indicator(t,:)==2),2)/task_num;
        drop_rate = size(find(indicator(t,:)==3),2)/task_num;
    end

    mode_num(t,:) = [local_rate,offloading_rate,drop_rate];
    
    disp(['在第',num2str(t),'轮:']);
    disp(['本地执行的移动设备占比: ',num2str(local_rate)]);
    disp(['卸载执行的移动设备占比: ',num2str(offloading_rate)]);
    disp(['任务丢弃的移动设备占比: ',num2str(drop_rate)]);

    % 移动设备电量迭代
    B(t+1,:) = B(t,:) - E(t,:) + e(t,:);
    % 时间片迭代
    t = t + 1;
end

%% 结果总结
disp('--------------迭代结束--------------');
disp(['本地执行的平均移动设备占比: ', num2str(mean(mode_num(:,1)))]);
disp(['卸载执行的平均移动设备占比: ', num2str(mean(mode_num(:,2)))]);
disp(['任务丢弃的平均移动设备占比: ', num2str(mean(mode_num(:,3)))]);