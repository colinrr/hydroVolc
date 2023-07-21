function beta_folch(
    Ri, ws
    )
    beta = 0.34 * (sqrt(2 * abs(Ri)) * ws) ^ -0.125
    if beta < 0.1
        beta = 0.1
    end
    if beta > 1
        beta = 1
    end  
    return beta
    end
