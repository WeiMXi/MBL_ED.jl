function lsr(eigval_sol)
    eigval_sol = real.(eigval_sol)
    sort!(eigval_sol)
    gap = Array{Float64}(undef, length(eigval_sol) - 1)
    r = Array{Float64}(undef, length(eigval_sol) - 2)
    for i in 1:length(eigval_sol) - 1
        gap[i] = eigval_sol[i + 1] - eigval_sol[i]
    end 
    
    for i in 1:length(eigval_sol) - 2
        r[i] = min(gap[i + 1], gap[i])/max(gap[i + 1], gap[i])
    end
    # gap_of_E
    mean(r)
end

function EE(eigvec_sol, left_chain, right_chain, nhil::Int)
    
    # entroy = Array{Float64}(undef, nhil)
    # for i in 1:nhil
    #     rho = sparse(left_chain, right_chain, eigvec_sol[:, i])
    #     rho_left = rho * rho'
    #     entroy[i] = log(tr(rho_left^2)) / (1 - 2)
    # end
    # mean(entroy[ceil(Int,nhil * 3/8):ceil(Int, nhil * 5/8)])

    # entroy = Array{Float64}(undef, nhil)
    # for i in 1:nhil
    #     rho = sparse(left_chain, right_chain, eigvec_sol[:, i])
    #     rho_left = rho * rho'
    #     rho_left = Array(rho_left)
    #     entroy[i] = real(-tr(rho_left * log(rho_left)))
    # end
    # mean(entroy)

    # entroy = Array{Float64}(undef, nhil)
    # # test = []
    # # my_flag::Bool = true
    # for i in 1:nhil
    #     rho = sparse(left_chain, right_chain, eigvec_sol[:, i])
    #     rho_left = rho * rho'
    #     rho_left = Symmetric(Array(rho_left))
    #     # eig_v = abs.(real(eigvals(rho_left)))
    #     eig_v = real.(eigvals(rho_left))
    #     # eig_v = LinearAlgebra.LAPACK.syev!('N', 'U', rho_left)
    #     entroy[i] = -sum(eig_v .* log.(abs.(eig_v)))
    #     if isnan(entroy[i])
    #         # push!(test, eig_v)
    #         entroy[i] = 0.0
    #     end
    #     # Symmetric[rho_left]
    #     # my_flag = my_flag && (sum(rho_left .== Symmetric(rho_left)) == length(rho_left))
    # end
    # mean(entroy[ceil(Int,nhil * 3/8):ceil(Int, nhil * 5/8)])    # 根据文献只取中间1/4
    # # test
    # # return my_flag

    entroy = Array{Float64}(undef, nhil)
    # test = []
    # my_flag::Bool = true
    for i in 1:nhil
        rho = sparse(left_chain, right_chain, eigvec_sol[:, i])
        

        #################################
        # GPU
        # rho = CuArray(Array(rho))
        # rho_left = rho * rho'
        # rho = nothing
        # # eig_v = CUSOLVER.syevd!( 'N', 'U', CuArray(Symmetric(Array(rho_left))) )
        # eig_v = CUSOLVER.syevd!( 'N', 'U', rho_left )
        # eig_v = Array(eig_v)
        #################################

        #################################
        # CPU
        rho_left = rho * rho'
        rho = nothing
        eig_v = eigvals( Symmetric(Array(rho_left)) )
        # eig_v = abs.(real(eigvals(rho_left)))
        # eig_v = real.(eigvals(rho_left))
        # eig_v = LinearAlgebra.LAPACK.syev!('N', 'U', rho_left)
        #################################

        entroy[i] = -sum(eig_v .* log.(abs.(eig_v)))
        if isnan(entroy[i])
            # push!(test, eig_v)
            entroy[i] = 0.0
        end
        # Symmetric[rho_left]
        # my_flag = my_flag && (sum(rho_left .== Symmetric(rho_left)) == length(rho_left))
    end
    mean(entroy[ceil(Int,nhil * 3/8):ceil(Int, nhil * 5/8)])    # 根据文献只取中间1/4
    # test
    # return my_flag
end