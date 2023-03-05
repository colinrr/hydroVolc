"""
alpha_cara docs:
    Ri: 
    z:
    D:
    Lm:
    R:
"""
function alpha_cara(
    Ri,z,D,Lm,R
    )
    if z/D < 15
        ajet = 1.1+0.004639652666*(z/D)^2-0.00019946078370*(z/D)^3
        aplume = 1.33-0.00343428985*(z/D)^2+0.00020643870*(z/D)^3
        dajet = (2*0.0046*(z/D)-3*0.0002*(z/D)^2)/D
        daplume = (0.0034*(z/D)*2+3*0.0002*(z/D)^2)/D
    
    else
        ajet = 2.45-1.05*exp(-0.00465*z/D)
        aplume = 1.42-4.42*exp(-0.2188*z/D)
        dajet = 1.05*(0.00465/D)*exp(-0.00465*z/D)
        daplume = 4.42*(0.2188/D)*exp(-0.2188*z/D)
        
    end
    atot = atot = ajet + ((aplume - ajet) * (z / Lm - 1)) / 4
    datot = dajet + ((daplume - dajet) * (z / Lm - 1)) / 4 + ((aplume - ajet) * (1 / Lm)) / 4
    
    
    
    alpha = 0.135 / 2 + (1 - 1 / atot) * Ri + (R * (datot / atot)) / 2
    
    if alpha < 0.05
        alpha = 0.05
    end
    if alpha > 0.17
        alpha = 0.17
    end
    
    return alpha
    end
