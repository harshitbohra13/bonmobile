def get_rho(h):
    T = 288.15 - 6.5*h/1000
    p = 101325*(1-0.0065*h/288.15)**5.2561
    rho = p/(287*T)
    return rho

print(get_rho(450))
