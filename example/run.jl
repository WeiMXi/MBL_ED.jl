using Distributed
using Dates
using Statistics
using JLD

num_of_cores_for_cal = 5
if length(workers()) < num_of_cores_for_cal
    addprocs(num_of_cores_for_cal - length(workers()) + 1)
end

@everywhere push!(LOAD_PATH,".")
@everywhere using MBL_ED

print("开始计算：")
ts = now()
println(ts)


the_ws = [i for i in 0.5:0.5:2.0]
append!(the_ws, [i for i in 2.25:0.25:5.0])
append!(the_ws, [5.5, 6.0, 7.0, 8.0, 9.0, 10.0])

r = give_empty_dataframe_for_format()
## 计算部分
# for n in [10]
for n in [10, 12]
    append!(r, primitive_run(tol_n = n, up = div(n, 2), W = the_ws, samples = 10))
end
# k = r[1]
# for i in 2:length(r)
#     append!(k, r[i])
# end
# save("/data/aybai/Dr._Li/fin_MBL_CUDA.jld", "result", k);


print("计算完成：")
tf = now()
println(tf)
println(Dates.value(tf - ts)/1000, "seconds used")

JLD.save("MBL_ED_test.jld", "result", r);
