using JLD   # 数据存取
using Gadfly    # 绘图
using DataFrames # 数据分类
using Cairo
using Statistics

# 需要计算置信区间的函数
function error_95(data)
    1.96*std(data)/sqrt(length(data))
end

function isequ(a, b)
    if a == b
        return true
    else
        return false
    end
end

r = load("MBL_ED_test.jld")
rDF = r["result"]
# rDF = filter(:E => E -> !(ismissing(E) || isnothing(E) || isnan(E)), rDF)
rDF = filter(( [:W, :Class] => (W, Class) -> (!( isequ(W, 8.5) || isequ(W, 9.0)|| isequ(W, 10.0)) || !( isequ(Class, "QP"))) ), rDF)

#########################################################
# 读取数据并封装至Dataframe
# r = load("fin_MBL_CPU.jld")
# d_c = Vector{String}(undef, length(r["result"]))
# # d_n = Vector{Int}(undef, length(r["result"]))
# d_n = Vector{String}(undef, length(r["result"]))
# d_w = Vector{Float64}(undef, length(r["result"]))
# d_r = Vector{Float64}(undef, length(r["result"]))
# d_e = Vector{Float64}(undef, length(r["result"]))
# for i in 1:length(r["result"])
#     d_c[i] = r["result"][i][1]
#     # d_n[i] = r["result"][i][2]
#     d_n[i] = string(r["result"][i][2])
#     d_w[i] = r["result"][i][3]
#     d_r[i] = r["result"][i][4]
#     d_e[i] = r["result"][i][5]
# end
# rDF = DataFrame(Class = d_c, N_tol = d_n, W = d_w, r = d_r, E = d_e)
# d_c, d_n, d_w, d_r, d_e = [nothing for i in 1:5]
# rDF = filter(:E => E -> !(ismissing(E) || isnothing(E) || isnan(E)), rDF)


# r = load("test_MBL_3.jld")
# d_c = Vector{String}(undef, length(r["result"]))
# # d_n = Vector{Int}(undef, length(r["result"]))
# d_n = Vector{String}(undef, length(r["result"]))
# d_w = Vector{Float64}(undef, length(r["result"]))
# d_r = Vector{Float64}(undef, length(r["result"]))
# d_e = Vector{Float64}(undef, length(r["result"]))
# for i in 1:length(r["result"])
#     d_c[i] = r["result"][i][1]
#     # d_n[i] = r["result"][i][2]
#     d_n[i] = string(r["result"][i][2])
#     d_w[i] = r["result"][i][3]
#     d_r[i] = r["result"][i][4]
#     d_e[i] = r["result"][i][5]
# end
# rDF_2 = DataFrame(Class = d_c, N_tol = d_n, W = d_w, r = d_r, E = d_e)
# d_c, d_n, d_w, d_r, d_e = [nothing for i in 1:5]
# rDF_2 = filter(:E => E -> !(ismissing(E) || isnothing(E) || isnan(E)), rDF_2)
# append!(rDF, rDF_2)
# rDF = filter(( [:W, :Class] => (W, Class) -> (!( isequ(W, 8.5) || isequ(W, 9.0)|| isequ(W, 10.0)) || !( isequ(Class, "QP"))) ), rDF)

#########################################################


# 数据分类，并计算平均值与置信区间
gdf = groupby(rDF, [:Class, :W, :N_tol])
fin_DF = combine(gdf, :E => error_95 , :E => mean, :r => error_95 , :r => mean)
fin_DF = groupby(fin_DF, :Class)

# 绘图
# 准备刻度
r_ticks = [0.38, 0.43, 0.48, 0.53]
e_ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
mysmoothing = 0.65 # 平滑度
# R Model plot
my_buffer = fin_DF[(Class = "R",)]
myp_Rr = plot(my_buffer, x=:W, y=:r_mean, color=:N_tol, shape=:N_tol
    , ymin = my_buffer.r_mean - my_buffer.r_error_95
    , ymax = my_buffer.r_mean + my_buffer.r_error_95
    , Geom.errorbar
    , Geom.point
    # , Geom.line
    , Guide.xticks(ticks=0:2:10)
    , Guide.yticks(ticks=r_ticks, label=true)
    , yintercept = [0.39, 0.53]
    , Geom.hline(style = [:dash, :dash], color=["gray", "gray"])
    , Guide.xlabel("W")
    , Guide.ylabel("Mean of r")
    , Guide.title("R Model")
    , layer( Stat.smooth(method=:loess, levels=[0.95], smoothing = mysmoothing), Geom.line )
    )

myp_RE = plot(my_buffer, x=:W, y=:E_mean, color=:N_tol, shape=:N_tol
    , ymin = my_buffer.E_mean - my_buffer.E_error_95
    , ymax = my_buffer.E_mean + my_buffer.E_error_95
    , Geom.errorbar
    , Geom.point
    # , Geom.line
    , Guide.xticks(ticks=0:2:10)
    , Guide.yticks(ticks=e_ticks)
    , Guide.xlabel("W")
    , Guide.ylabel("Entroy")
    , Guide.title("R Model")
    , layer( Stat.smooth(method=:loess, levels=[0.95], smoothing = mysmoothing), Geom.line )
    )

# img = SVG("r_w_test.svg", 6inch, 4inch)
# draw(img, myp)

# QP Model plot
my_buffer = fin_DF[(Class = "QP",)]

myp_QPr = plot(my_buffer, x=:W, y=:r_mean, color=:N_tol, shape=:N_tol
    , ymin = my_buffer.r_mean - my_buffer.r_error_95
    , ymax = my_buffer.r_mean + my_buffer.r_error_95
    , Geom.errorbar
    , Geom.point
    # , Geom.line
    , yintercept = [0.39, 0.53]
    , Geom.hline(style = [:dash, :dash], color=["gray", "gray"])
    , Guide.xticks(ticks=0:2:8)
    , Guide.yticks(ticks=r_ticks)
    , Guide.xlabel("W")
    , Guide.ylabel("Mean of r")
    , Guide.title("QP Model")
    , layer( Stat.smooth(method=:loess, levels=[0.95], smoothing = mysmoothing), Geom.line )
    )

myp_QPE = plot(my_buffer, x=:W, y=:E_mean, color=:N_tol, shape=:N_tol
    , ymin = my_buffer.E_mean - my_buffer.E_error_95
    , ymax = my_buffer.E_mean + my_buffer.E_error_95
    , Geom.errorbar
    , Geom.point
    # , Geom.line
    , Guide.xticks(ticks=0:2:8)
    , Guide.yticks(ticks=e_ticks)
    , Guide.xlabel("W")
    , Guide.ylabel("Entroy")
    , Guide.title("QP Model")
    , layer( Stat.smooth(method=:loess, levels=[0.95], smoothing = mysmoothing), Geom.line )
    )

my_buffer = nothing

# img = SVG("E_w_test.svg", 6inch, 4inch)
# draw(img, myp)


fin_pl = gridstack([myp_QPE myp_RE; myp_QPr myp_Rr])

# img = SVG("5_result_CUDA_test.svg", 12inch, 8inch)
img = PDF("fin_MBL_CPU.pdf", 12inch, 8inch)
# img_PNG = PNG("result_test.svg", 12inch, 8inch)
draw(img, fin_pl)
