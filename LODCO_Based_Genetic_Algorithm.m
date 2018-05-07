clc, clear
%% ================= Simulation of LODCO-Based Genetic Algorithm =================

%% 基本参数设置
k = 1e-28;                        % 有效开关电容
tau = 0.002;                      % 时间片长度(s)
tau_d = 0.002;                    % 计算任务执行时间的deadline(s)
phi = 0.002;                      % 任务丢弃的惩罚项权重(s)
little_phi = 0.002;               % 用小写的phi表示任务被卸载执行的奖励项
omega = 1e6;                      % 服务器带宽(Hz)
sigma = 1e-13;                    % 接收端的噪声功率(W)
p_tx_max = 1;                     % 移动设备最大传输功率(W)
f_max = 1.5e9;                    % 移动设备最大CPU时钟周期频率(Hz)
E_max = 0.002;                    % 电池允许的最大放电量(J)
L = 1000;                         % 一项计算任务的大小(bit)
X = 737.5;                        % 移动设备执行一项计算任务所需的时钟周期个数
W = 737500;                       % 移动设备本地执行一项计算任务所需的时钟周期个数(L*X)
g0 = 1e-4;                        % 路径损失常数(dB转化之后的数值比)
d0 = 1;                           % 服务器和移动设备之间的相对距离(m)

%% 变量控制
N = 10;                           % 移动设备数目
M = 5;                            % MEC服务器个数
T = 150;                          % 时间片个数
E_min = 0.02e-3;                  % 能量使用下界(J)
V = 1e-5;                         % LODCO中penalty项的权重(J^2/second)
rho = 0.7;                        % 计算任务抵达的概率
E_H_max = 48e-6;                  % 收集的能量服从的均匀分布上限(J)
eps = 0.1;                        % 用于贪心算法的决策

% 实际能耗上限
E_max_hat = min(max(k*W*(f_max^2), p_tx_max*tau), E_max);
theta = E_max_hat + V*phi/E_min;        % 扰动参数

%% 中间变量存储           
B = zeros(T, N);                        % N个移动设备的实际电量
B_hat = zeros(T, N);                    % N个移动设备的虚拟电量
e = zeros(T, N);                        % N个移动设备的能量收集
indicator = zeros(T, N);                % 用1,2,3,4分别表示本地执行、卸载执行、drop(前三个均意味着有任务到达)以及没有任务到达
f = zeros(T, N);                        % N个移动设备本地执行的频率(之后不会用到)
p = zeros(T, N);                        % N个移动设备卸载执行的传输功率(之后不会用到)
local_execution_delay = zeros(T, N);    % N个移动设备的local execution delay
remote_execution_delay = zeros(T, N);   % N个移动设备的remote execution delay
cost = zeros(T, N);                     % N个移动设备的execution cost(最后确定)
E_local = zeros(T, N);                  % N个移动设备本地执行的能耗
E_remote = zeros(T, N);                 % N个移动设备卸载执行的能耗
E_all = zeros(T, N);                    % N个移动设备的总能耗
mode_num = zeros(T,3);                  % 每一列分别表示每轮中本地执行、远程执行及任务丢弃的比率(分母为N减去没有任务到达的)

% 关闭fslove的输出
opt = optimset('Display', 'off');

t = 1;
while t <= T
    % 用一个不定长的、3列的矩阵描述"移动设备i-MEC服务器j-i到j的最小传输延迟"
    map = [];
    % 存储每一个MEC服务器连接到的移动设备个数
    flags = zeros(M, 1);
    % 分别保存当前移动设备到各MEC服务器的J_s值(就是延迟)
    J_s_matrix = zeros(N, M);
    % 分别保存当前移动设备到各MEC服务器的能耗
    E_remote_matrix = zeros(N, M);
    % 分别保存当前移动设备到各MEC服务器的最佳传输功率
    p_matrix = zeros(N, M);
    % 分别存储每一个移动设备的本地执行的延迟J_m供二次决策时使用
    J_m = zeros(N, 1);
    % 产生虚拟电量队列值
    B_hat(t,:) = B(t,:) - theta;

    % 假设一般情况下不需要借助键值对(即直接根据intlinporg求解)
    useKeyValuePair = 0;
    for i = 1:N
        %% 对每一个移动设备，阶段初始化
        % 以伯努利分布产生计算任务
        zeta = binornd(1, rho);
        if zeta == 0
            % 没有计算任务产生
            indicator(t, i) = 4;
            f(t, i) = 0; p(t, i) = 0;
        else
            %% 求解optimal energy harvesting e*
            % 产生E_H_t
            E_H_t = unifrnd(0, E_H_max);
            if B_hat(t,i) <= 0                    % 初始值为0，因此无需讨论大于0的情形
                e(t, i) = E_H_t;
            end

            %% 求解P_ME
            f_L = max(sqrt(E_min/(k*W)), W/tau_d);
            f_U = min(sqrt(E_min/(k*W)), f_max);
            if f_L <= f_U
                % P_ME有解
                f0 = power(V/(-1*B_hat(t,i)*k), 1/3);
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
                    disp(['P_ME电量不足![此时t为', num2str(t), ']']);
                    J_m(i) = inf;    % 设置为inf可以保证一定比phi小
                    useKeyValuePair = 1;
                else
                    % 计算此时的J_m(只考虑延迟，不计算能耗)
                    J_m(i) = W/f(t, i);
                end

            else
                disp(['P_ME无解![此时t为', num2str(t), ']']);
                % 因此indicator(t, 1) = 0不变
                J_m(i) = inf;
                useKeyValuePair = 1;
            end

            %% 求解P_SE
            % 随机产生服务器和移动设备的距离(限定在0 ~ 60之内)
            D = unifrnd(0, 70, N, M);
            % 服从lambda=1的指数分布的小尺度衰落信道功率收益
            gamma = exprnd(1, N, M);
            % 从任意移动设备到任意服务器的信道功率增益
            h = g0*gamma.*power(d0./D, 4);

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
                    tmp = fsolve(y, [0.001, 1], opt);
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
                    p0 = real(max(fsolve(y, [0.001, 1], opt)));
                    if p_U < p0
                        p_matrix(i, j) = p_U;
                    elseif p_L > p0 && B_hat(t,i) < 0
                        p_matrix(i, j) = p_L;
                    elseif p_L <= p0 && p_U >= p0 && B_hat(t,i) < 0
                        p_matrix(i, j) = p0;
                    end
                    % 计算achievable rate
                    r = calAchieveRate(tmp_h, p_matrix(i, j), omega, sigma);
                    % 计算此时的能耗
                    E_remote(t, i) = p_matrix(i, j) * L/r;
                    E_remote_matrix(i, j) = E_remote(t, i);
                    if E_remote(t, i) >= B(t, i)
                        disp(['P_SE电量不足![此时t为', num2str(t), ',移动设备编号为', num2str(i), ',MEC服务器编号为', num2str(j), '].']);
                        J_s = inf;
                        useKeyValuePair = 1;
                    else
                        % 计算此时的J_s(只计算延迟，不考虑能耗)
                        J_s = L/r;
                    end
                else
                    disp(['P_SE无解![此时t为', num2str(t), ',移动设备编号为', num2str(i), '].']);
                    J_s = inf;
                    useKeyValuePair = 1;
                end
                J_s_matrix(i,j) = J_s;
            end
            % 计算此时的execution delay
            remote_execution_delay(t, i) = min(J_s_matrix(i,:));
            % 保存最佳的execution delay及其对应的服务器编号
            [J_s_best, j_best] = min(J_s_matrix(i,:));
            E_remote(t, i) = E_remote_matrix(i, j_best);
            %% 为第i个移动设备选取最佳模式
            [~, mode] = min([J_m(i), J_s_best, phi]);
            indicator(t, i) = mode;
            if mode == 2
                map = [map;[i,j_best,J_s_best]];
            end
        end
    end


    %% 求解过程
    if useKeyValuePair == 0
        % 亦采用intlinprog求解        
        % 定义目标函数f
        goal = zeros(1,N*(M+2));
        for i = 1:10
            goal(7*i-6:7*i-5) = [local_execution_delay(t,i)-B_hat(t,i)*E_local(t,i),phi];
            goal(7*i-4:7*i) = J_s_matrix(i,:)-B_hat(t,i)*E_remote_matrix(i,:)-little_phi;
        end
        fitnessGoal = @(x)-goal*[x(1),x(2),x(3),x(4),x(5),x(6),x(7),...
        x(8),x(9),x(10),x(11),x(12),x(13),x(14),...
        x(15),x(16),x(17),x(18),x(19),x(20),x(21),...
        x(22),x(23),x(24),x(25),x(26),x(27),x(28),...
        x(29),x(30),x(31),x(32),x(33),x(34),x(35),...
        x(36),x(37),x(38),x(39),x(40),x(41),x(42),...
        x(43),x(44),x(45),x(46),x(47),x(48),x(49),...
        x(50),x(51),x(52),x(53),x(54),x(55),x(56),...
        x(57),x(58),x(59),x(60),x(61),x(62),x(63),...
        x(64),x(65),x(66),x(67),x(68),x(69),x(70)]';
        % 定义nvars
        nvars = N*(M+2);
        % 定义A, b, lb, ub
        A = [0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0;
             0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0;
             0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0;
             0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0;
             0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1;
             1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0;
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1];
        b = [4;4;4;4;4;1;1;1;1;1;1;1;1;1;1];
        lb = zeros(N*(M+2),1);
        ub = ones(N*(M+2),1);
        % 返回计算结果 (system operation)
        options = gaoptimset('ParetoFraction',0.3,'PopulationSize',100,'Generations',200,'StallGenLimit',200,'TolFun',1e-100,'PlotFcns',@gaplotpareto);
        so = gamultiobj(fitnessGoal,nvars,A,b,[],[],lb,ub,options);
        for i = 1:10
            pos = find(so(7*i-6:7*i)==1);
            if pos == 1
                indicator(t,i) = 1;
            elseif pos == 2
                indicator(t,i) = 3;
            else
                indicator(t,i) = 2;
            end
        end
    else
        % 采用键值对求解
        %% 为那些选择卸载执行的任务分配服务器(或另选模式)
        % UB = f_server_max*tau_d/(L*X),取f_server_max=f_max则约为4.0678
        UB = 4;
        while ~isempty(map)
            % 找到拥有最小传输延迟的移动设备-MEC服务器对
            [min_Js,index] = min(map(:,3));
            % 找到最小延迟对应的min_i和min_j
            min_i = map(index,1);
            min_j = map(index,2);
            
            % 此时只考虑卸载执行，即不比较卸载执行是否为三者最佳
            if rand() <= eps
                if flags(min_j) <= UB
                    % 从map中删除该键值对并同步一系列共同维护的变量
                    map(index,:) = [];
                    % 对应的MEC服务器自增1(该操作的位置发生了变化，论文此处需要修改)
                    flags(min_j) = flags(min_j) + 1;
                    % 将J_s_matrix(min_i,min_j)设为inf
                    J_s_matrix(min_i,min_j) = inf;
                else
                    if min(J_s_matrix(min_i,:)) ~= inf
                        % 当该移动设备还有能选的服务器的时候，不断找能找到的最小的那个
                        % 此处论文需要补充！找到最小的那个之后，覆盖map中该移动设备选取的服务器，并覆盖对应的Js最小值
                        [min_Js_second, min_j_second] = min(J_s_matrix(min_i,:));
                        map(index,2:3) = [min_j_second,min_Js_second];
                        % 返回最外层的while，重新开始找最小的Js
                        continue;
                    else
                        % 没有服务器可以选了，只能在另外两种mode中选取
                        % 重新设置指示变量
                        [~, mode] = min([J_m(min_i), inf, phi]);
                        indicator(t, i) = mode;
                        % 从map中删除该键值对并同步一系列共同维护的变量(此处论文需要补充)
                        map(index,:) = [];
                        %{
                        ==min_i已绝无再出现的可能，因此没有必要再将J_s_matrix(min_i,min_j)设为inf==
                        J_s_matrix(min_i,min_j) = inf;
                        %}
                    end
                end
            else
                [~, mode] = min([J_m(i), J_s_best, phi]);
                if mode == 2
                    % 当前最优模式仍为卸载执行
                    if flags(min_j) <= UB
                        % 从map中删除该键值对并同步一系列共同维护的变量
                        map(index,:) = [];
                        % 对应的MEC服务器自增1(该操作的位置发生了变化，论文此处需要修改)
                        flags(min_j) = flags(min_j) + 1;
                        % 将J_s_matrix(min_i,min_j)设为inf
                        J_s_matrix(min_i,min_j) = inf;
                    else
                        if min(J_s_matrix(min_i,:)) ~= inf
                            % 当该移动设备还有能选的服务器的时候，不断找能找到的最小的那个
                            % 此处论文需要补充！找到最小的那个之后，覆盖map中该移动设备选取的服务器，并覆盖对应的Js最小值
                            [min_Js_second, min_j_second] = min(J_s_matrix(min_i,:));
                            map(index,2:3) = [min_j_second,min_Js_second];
                            % 返回最外层的while，重新开始找最小的Js
                            continue;
                        else
                            % 没有服务器可以选了，只能在另外两种mode中选取
                            % 重新设置指示变量
                            [~, mode] = min([J_m(min_i), inf, phi]);
                            indicator(t,i) = mode;
                            % 从map中删除该键值对并同步一系列共同维护的变量(此处论文需要补充)
                            map(index,:) = [];
                            %{
                            ==min_i已绝无再出现的可能，因此没有必要再将J_s_matrix(min_i,min_j)设为inf==
                            J_s_matrix(min_i,min_j) = inf;
                            %}
                        end
                    end
                else
                    % 设置新的最优模式，并从map中删除本键值对及维护的变量
                    indicator(t,i) = mode;
                    map(index,:) = [];
                end
            end
        end
    end

    % 计算每一个移动设备的execution cost
    cost(t,indicator(t,:)==1) = local_execution_delay(t,indicator(t,:)==1);
    cost(t,indicator(t,:)==2) = remote_execution_delay(t,indicator(t,:)==2);
    cost(t,indicator(t,:)==3) = phi;
    % 计算每一个移动设备的总能耗
    E_all(t,indicator(t,:)==1) = E_local(t,indicator(t,:)==1);
    E_all(t,indicator(t,:)==2) = E_remote(t,indicator(t,:)==2);
    E_all(t,indicator(t,:)==3) = 0;

    % 中间结果输出
    task_num = N - size(find(indicator(t,:)==4),2);
    local_rate = size(find(indicator(t,:)==1),2)/task_num;
    offloading_rate = size(find(indicator(t,:)==2),2)/task_num;
    drop_rate = size(find(indicator(t,:)==3),2)/task_num;

    mode_num(t,:) = [local_rate,offloading_rate,drop_rate];
    
    disp(['在第',num2str(t),'轮:']);
    disp(['本地执行的移动设备占比: ',num2str(local_rate)]);
    disp(['卸载执行的移动设备占比: ',num2str(offloading_rate)]);
    disp(['任务丢弃的移动设备占比: ',num2str(drop_rate)]);
    disp('-----------------------------------');

    % 移动设备电量迭代
    B(t+1,:) = B(t,:) - E_all(t,:) + e(t,:);
    % 时间片迭代
    t = t + 1;
end

%% 结果总结
disp('--------------迭代结束--------------');
disp(['本地执行的平均移动设备占比: ', num2str(mean(mode_num(:,1)))]);
disp(['卸载执行的平均移动设备占比: ', num2str(mean(mode_num(:,2)))]);
disp(['任务丢弃的平均移动设备占比: ', num2str(mean(mode_num(:,3)))]);