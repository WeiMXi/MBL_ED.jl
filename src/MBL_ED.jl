# baiaiyu
# 2021,Jul,15
# Many-Body Localization Transition https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.075702
# programe for application
# 请君洗耳 听我之言

# bitstring() 获取二进制结构
# 稀疏矩阵
# ↓ 实对称矩阵在julia
# https://discourse.julialang.org/t/are-eigenvectors-of-real-symmetric-or-hermitian-matrix-computed-from-eigen-all-orthogonal/60893/12
# ↓ LAPACK在julia
# LinearAlgebra.LAPACK.syevr!

module MBL_ED
    
export cal_one_task, primitive_run, give_tasks, auto_manger, give_empty_dataframe_for_format

using SparseArrays  # 稀疏矩阵支持
using LinearAlgebra # QR法得特征值和特征向量
# using KrylovKit     # lanczos 直接得到特征值和特征向量
# using Arpack
using Statistics    # 平均值与标准差
# using Gadfly        # 绘图
using Distributed   #   并行
using JLD   # 保存数据
using DataFrames    # 数据整合
using Dates # 时间

######### GPU ##########
# using CUDA.CUSOLVER
# using CUDA
########################

include("./my_hamiltonion.jl")
include("./phy_quanlities.jl")
include("./mymath.jl")
include("./parallel_manager.jl")    # 并行计算管理器
include("./mydata.jl")    # 数据整合


# 封装计算任务
struct QPTask
    H::Matrix{Float64}
    states::Vector{Int}
    left_chain::Vector{Int}
    right_chain::Vector{Int}
    nhil::Int
    n_tol::Int
    W::Float64
    omega::Float64
    phi::Float64
end

struct RTask
    H::Matrix{Float64}
    states::Vector{Int}
    left_chain::Vector{Int}
    right_chain::Vector{Int}
    nhil::Int
    n_tol::Int
    W::Float64
    omega::Float64
end

function cal_one_task(one_R_task::RTask)
    St = 0.5*(one_R_task.n_tol*log(2) - 1)
    # HH = Array{Float64}(undef, one_R_task.nhil, one_R_task.nhil)
    # HH .= one_R_task.H
    HH = copy(one_R_task.H)
    wave = fill_Sz(one_R_task.states, one_R_task.n_tol, one_R_task.nhil, one_R_task.W, one_R_task.omega)
    for i in 1:one_R_task.nhil
        @inbounds HH[i, i] = HH[i, i] + wave[i]
    end

    empty_df = give_empty_dataframe_for_format()

    #################################
    # CPU
    eig_sol = eigen(Symmetric(HH))
    push!( empty_df, ["R", string(one_R_task.n_tol), one_R_task.W, lsr(eig_sol.values), EE(eig_sol.vectors, one_R_task.left_chain, one_R_task.right_chain, one_R_task.nhil)/St] )
    #################################

    #################################
    # GPU
    # eig_sol = CUSOLVER.syevd!( 'V', 'U', CuArray(Symmetric(HH)) )
    # #"R", one_R_task.n_tol, one_R_task.W, lsr(Array(eig_sol[1])), EE(Array(eig_sol[2]), one_R_task.left_chain, one_R_task.right_chain, one_R_task.nhil)/St
    # push!( empty_df, ["R", string(one_R_task.n_tol), one_R_task.W, lsr(eig_sol.values), EE(eig_sol.vectors, one_R_task.left_chain, one_R_task.right_chain, one_R_task.nhil)/St] )
    #################################
end

function cal_one_task(one_QP_task::QPTask)
    St = 0.5*(one_QP_task.n_tol*log(2) - 1)
    # HH = Array{Float64}(undef, one_QP_task.nhil, one_QP_task.nhil)
    # HH .= one_QP_task.H
    HH = copy(one_QP_task.H)
    wave = fill_Sz(one_QP_task.states, one_QP_task.n_tol, one_QP_task.nhil, one_QP_task.W, one_QP_task.omega, one_QP_task.phi)
    for i in 1:one_QP_task.nhil
        @inbounds HH[i, i] = HH[i, i] + wave[i]
    end

    empty_df = give_empty_dataframe_for_format()

    #################################
    # CPU
    eig_sol = eigen(Symmetric(HH))
    push!(empty_df, ["QP", string(one_QP_task.n_tol), one_QP_task.W, lsr(eig_sol.values), EE(eig_sol.vectors, one_QP_task.left_chain, one_QP_task.right_chain, one_QP_task.nhil)/St] )
    #################################


    #################################
    # GPU
    # eig_sol = CUSOLVER.syevd!( 'V', 'U', CuArray(Symmetric(HH)) )
    # #"QP", one_QP_task.n_tol, one_QP_task.W, lsr(Array(eig_sol[1])), EE(Array(eig_sol[2]), one_QP_task.left_chain, one_QP_task.right_chain, one_QP_task.nhil)/St
    # push!(empty_df, ["QP", one_QP_task.n_tol, one_QP_task.W, lsr(eig_sol.values), EE(eig_sol.vectors, one_QP_task.left_chain, one_QP_task.right_chain, one_QP_task.nhil)/St] )
    #################################
end

function primitive_run(
            ; tol_n::Int = 10
            , up::Int = 5
            , J = 1.0
            , Jp = 1.0
            , Jz = 1.0
            , W = 0.5:0.5:8.0
            , k = (sqrt(5) - 1)/2
            , samples::Int = 10
            )

    omega = 2*pi*k
    # Hilbert space 维数
    nhil::Int = div(factorial(big(tol_n)), (factorial(big(up))*factorial(big(tol_n - up))))
    states = Array{Int}(undef, nhil)
    # 得到所有要求的状态
    all_required_states!(states, nhil, tol_n, up)
    # 得到左半链和右半链的状态
    left_states = Array{Int}(undef, nhil)
    right_states = Array{Int}(undef, nhil)
    for i in 1:nhil
        left_states[i], right_states[i] = left_right_half_chain(states[i], tol_n)
    end
    # 得到不含随机项的哈密顿量
    H = fill_H(nhil, states, tol_n, J, Jp, Jz)


    cal_tasks_array = []
    for w in W
        for i in 1:samples
            push!(cal_tasks_array, RTask(H, states, left_states, right_states, nhil, tol_n, w, omega))
            push!(cal_tasks_array, QPTask(H, states, left_states, right_states, nhil, tol_n, w, omega, (2.0*pi*(rand() - 0.5))))
        end
    end
    
    # println("getalltasks")
    auto_manger(cal_tasks_array)

    # one_R_task = RTask(H, states, left_states, right_states, nhil, tol_n, W, omega)
    # cal_one_task(cal_tasks_array[2])
    # cal_tasks_array
end

"""

    give_tasks(
    ; tol_n::Int = 10
    , up::Int = 5
    , J = 1.0
    , Jp = 1.0
    , Jz = 1.0
    , W = 0.5:0.5:8.0
    , k = (sqrt(5) - 1)/2
    , samples::Int = 10
    )

DOI: 10.1103/PhysRevLett.119.075702  
给出如下哈密顿量对应的计算任务, 一半为R模型一半为QP模型, 非周期性边界条件

``
\\begin{aligned}
H^{\\mathrm{QP} / R}=& J \\sum_{i=1}^{L-1}\\left(S_{i}^{x} S_{i+1}^{x}+S_{i}^{y} S_{i+1}^{y}\\right)+J_{z} \\sum_{i=1}^{L-1} S_{i}^{z} S_{i+1}^{z} \\\\
&+\\sum_{i=1}^{L} W \\cos \\left(2 \\pi k i+\\phi_{i}^{\\mathrm{QP} / R}\\right) S_{i}^{z} \\\\
&+J^{\\prime} \\sum_{i=1}^{L-2}\\left(S_{i}^{x} S_{i+2}^{x}+S_{i}^{y} S_{i+2}^{y}\\right)
\\end{aligned}
``
# Arguments
- `tol_n::Int`: 一维链长度
- `up::Int`: Sz向上的数量
- `samples::Int`: 对于每组参数的采样数

# Examples

```julia-repl
julia> tasks = give_tasks();
julia> tasks[1]
MBL_ED.RTask([2.75 0.5 … 0.0 0.0; 0.5 1.75 … 0.0 0.0; … ; 0.0 0.0 … 1.75 0.5; 0.0 0.0 … 0.5 2.75], [31, 47, 55, 59, 61, 62, 79, 87, 91, 93  …  15361, 15362, 15364, 15368, 15376, 15392, 15424, 15488, 15616, 15872], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  121, 121, 121, 121, 121, 121, 121, 122, 123, 125], [32, 48, 56, 60, 62, 63, 80, 88, 92, 94  …  2, 3, 5, 9, 17, 33, 65, 1, 1, 1], 2002, 14, 0.5, 3.883222077450933)

julia> tasks[2]
MBL_ED.QPTask([2.75 0.5 … 0.0 0.0; 0.5 1.75 … 0.0 0.0; … ; 0.0 0.0 … 1.75 0.5; 0.0 0.0 … 0.5 2.75], [31, 47, 55, 59, 61, 62, 79, 87, 91, 93  …  15361, 15362, 15364, 15368, 15376, 15392, 15424, 15488, 15616, 15872], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  121, 121, 121, 121, 121, 121, 121, 122, 123, 125], [32, 48, 56, 60, 62, 63, 80, 88, 92, 94  …  2, 3, 5, 9, 17, 33, 65, 1, 1, 1], 2002, 14, 0.5, 3.883222077450933, 2.0325482732260935)
```
"""
function give_tasks(
        ; tol_n::Int = 10
        , up::Int = 5
        , J = 1.0
        , Jp = 1.0
        , Jz = 1.0
        , W = 0.5:0.5:8.0
        , k = (sqrt(5) - 1)/2
        , samples::Int = 10
        )

    omega = 2*pi*k
    # Hilbert space 维数
    nhil::Int = div(factorial(big(tol_n)), (factorial(big(up))*factorial(big(tol_n - up))))
    states = Array{Int}(undef, nhil)
    # 得到所有要求的状态
    all_required_states!(states, nhil, tol_n, up)
    # 得到左半链和右半链的状态
    left_states = Array{Int}(undef, nhil)
    right_states = Array{Int}(undef, nhil)
    for i in 1:nhil
    left_states[i], right_states[i] = left_right_half_chain(states[i], tol_n)
    end
    # 得到不含随机项的哈密顿量
    H = fill_H(nhil, states, tol_n, J, Jp, Jz)


    cal_tasks_array = []
    for w in W
        for i in 1:samples
            push!(cal_tasks_array, RTask(H, states, left_states, right_states, nhil, tol_n, w, omega))
            push!(cal_tasks_array, QPTask(H, states, left_states, right_states, nhil, tol_n, w, omega, (2.0*pi*(rand() - 0.5))))
        end
    end

    # auto_manger(cal_tasks_array)

    # one_R_task = RTask(H, states, left_states, right_states, nhil, tol_n, W, omega)
    # cal_one_task(cal_tasks_array[2])
    cal_tasks_array
end

# 准备不含随机项的哈密顿矩阵等计算所需数据与模型参数
# 包括tol_n, up, nhil
function prepare_H(
    ; tol_n::Int = 10
    , up::Int = 5
    , J = 1.0
    , Jp = 1.0
    , Jz = 1.0
    )

    # Hilbert space 维数
    nhil::Int = div(factorial(big(tol_n)), (factorial(big(up))*factorial(big(tol_n - up))))
    states = Array{Int}(undef, nhil)
    # 得到所有要求的状态
    all_required_states!(states, nhil, tol_n, up)
    # 得到左半链和右半链的状态
    left_states = Array{Int}(undef, nhil)
    right_states = Array{Int}(undef, nhil)
    for i in 1:nhil
    left_states[i], right_states[i] = left_right_half_chain(states[i], tol_n)
    end
    # 得到不含随机项的哈密顿量
    H = fill_H(nhil, states, tol_n, J, Jp, Jz)

    file_name = "$tol_n" * "_$up"* "_.jld"
    save(file_name, "H", H, "states", states, "left_chain", left_states,"right_chain", right_states,"tol_n", tol_n, "up", up, "nhil", nhil)
    
    nothing
    # H
end

# 使用以保存的哈密顿量计算
struct spTask
    Hjld::String
    omega::Float64
    W::Float64
    phi::Float64
end

# 使用已经储存的哈密顿量计算
function cal_one_sptask(one_sptask::spTask)
    
    omega = 2*pi*k

    # 载入数据
    data = load(one_sptask.Hjld)

    cal_one_task(
        QPTask(data["H"], data["states"], data["left_chain"], data["right_chain"], data["nhil"], data["n_tol"], W, one_sptask.omega, (2.0*pi*(rand() - 0.5)))
    )


end

function run_sptasks(::Vector{String},)
    
end



end # module end