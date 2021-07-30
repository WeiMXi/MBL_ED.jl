#返回左右半链对应的序号，以Int的形式
function left_right_half_chain(the_state, tol_n)
    tr_bits = div(tol_n, 2)
    left_half::Int = (the_state >>> tr_bits)
    right_half::Int = (the_state - (left_half << tr_bits))
    # Julia数组从1开始
    left_half = left_half + 1
    right_half = right_half + 1
    left_half, right_half
end
