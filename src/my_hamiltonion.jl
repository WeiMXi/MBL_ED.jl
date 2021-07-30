# 计算哈密顿矩阵
function fill_H(nhil, states, tol_n, J, Jp, Jz)
    # 哈密顿量, 稀疏矩阵
    # H = spzeros(Float64, nhil, nhil)
    # 哈密顿量，full ED
    H = zeros(Float64, nhil, nhil)
    for i in 1:nhil
        ######################################################################
        # 计算1到tol_n - 2个格点的 S_-(i)S_+(j) + S_-(j)S_+(i) 值
        # 包括最近邻和次近邻
        for x in 1:tol_n - 2
            # 最近邻
            index = fill_up_down(states, i, x, x + 1)
            if index != 0.0
                H[i, index] = H[i, index] + 0.5*J
            end
            # 次近邻
            index = fill_up_down(states, i, x, x + 2)
            if index != 0.0
                H[i, index] = H[i, index] + 0.5*Jp
            end
        end
        # 补充计算tol_n - 1到tol_n格点的 S_-(i)S_+(j) + S_-(j)S_+(i) 值
        # 只有最近邻
        index = fill_up_down(states, i, tol_n - 1, tol_n)
        if index != 0.0
            H[i, index] = H[i, index] + 0.5*J
        end
        ######################################################################
        # 计算含Sz项的值，即对角线上的值（不含随机项
        H[i, i] = H[i, i] + fill_SzSz(states[i], tol_n, Jz)
        ######################################################################
    end
    H
end

# 给出态states[label] 的值: S_-(i)S_+(j) + S_-(j)S_+(i) 
# 为0则给出0
# 不为0则给出作用后态的label
function fill_up_down(states, label, i, j)
    the_state = states[label]
    i = i - 1 #从第一号平移i - 1位
    j = j - 1 #从第一号平移j - 1位
    i_state = btest(the_state, i)
    j_state = btest(the_state, j)
    # 态相同直接归零
    if i_state == j_state
        return 0.0
    end
    # 态不同则寻找作用后的态的label
    ## 先给出作用后的态
    if (i_state == 0 && j_state == 1)   # 作用后为i+ j-
        fin_state = the_state + 2^i - 2^j
    elseif (i_state == 1 && j_state == 0)   # 作用后为i- j+
        fin_state = the_state - 2^i + 2^j
    end
    searchsortedfirst(states, fin_state)
end

# 计算态 the_state 所有含Sz的项的值（不含随机项
# 返回计算结果
function fill_SzSz(the_state, tol_n, Jz)
    fin_result = 0.0
    for i in 1:tol_n - 1
        i_state = btest(the_state, i - 1)
        j_state = btest(the_state, i)
        #######################################
        # 计算Jz*Sz(i)*Sz(i + 1)
        if i_state == j_state
            fin_result = fin_result + 0.25*Jz
        else
            fin_result = fin_result - 0.25*Jz
        end
        #######################################
    end
    fin_result
end

# R 模型
# 给出含随机项的值
# 返回计算结果
function fill_Sz(states, tol_n, nhil::Int, W, omega)
    fin_result = zeros(Float64, nhil)
    rand_phi = rand(tol_n)
    for j in 1:nhil
        the_state::Int = states[j]
        for i in 1:tol_n
            i_state::Int = btest(the_state, i - 1)
        # 计算Sz(i)*w*cos(2pi*k*i + phi)
            if i_state == 1
                @inbounds fin_result[j] = fin_result[j] + 0.5*W*cos(omega*i + 2.0*pi*(rand_phi[i] - 0.5))
            else
                @inbounds fin_result[j] = fin_result[j] - 0.5*W*cos(omega*i + 2.0*pi*(rand_phi[i] - 0.5))
            end
        end
    end
    fin_result
end

# 给出含随机项的值,待乘正弦波因子
# 返回计算结果
# function fill_Sz(states, tol_n, nhil::Int)
#     fin_result = zeros(Float64, nhil)
#     for j in 1:nhil
#         the_state = states[j]
#         for i in 1:tol_n
#             i_state = btest(the_state, i - 1)
#         # 计算Sz(i)*w*cos(2pi*k*i + phi)
#             if i_state == 1
#                 fin_result[j] = fin_result[j] + 0.5
#             else
#                 fin_result[j] = fin_result[j] - 0.5
#             end
#         end
#     end
#     fin_result
# end

# QP 模型
# 给出含随机项的值
# 返回计算结果
function fill_Sz(states, tol_n, nhil::Int, W, omega, phi)
    fin_result = zeros(Float64, nhil)
    for j in 1:nhil
        the_state::Int = states[j]
        for i in 1:tol_n
            i_state::Int = btest(the_state, i - 1)
        # 计算Sz(i)*w*cos(2pi*k*i + phi)
            if i_state == 1
                @inbounds fin_result[j] = fin_result[j] + 0.5*W*cos(omega*i + phi)
            else
                @inbounds fin_result[j] = fin_result[j] - 0.5*W*cos(omega*i + phi)
            end
        end
    end
    fin_result
end

# 给出Sz为上数目为up的所有态，共有nhil个
function all_required_states!(states, nhil::Int, tol_n::Int, up::Int)
    num_alls::Int = 2^tol_n
    num_count = 0
    for i::Int in 1:num_alls
        num_up = 0
        for j in 0:tol_n-1
            if btest(i, j) == 1
                num_up = num_up + 1
            end
        end
        if num_up == up
            num_count = num_count + 1
            states[num_count] = i
            # println(bitstring(i))
        end
    end
    
    # 错误处理
    if num_count != nhil
        throw("dimension of Hilbert space mismatch nhil:$nhil num_count$num_count")
    end
end



# 返回s态第j位结果，上为1，下为0
function btest(s::Int, j::Int)
    u0::UInt = 1
    (s >>> j) & u0
end