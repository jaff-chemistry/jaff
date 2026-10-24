# ABOUTME: Template-engine regression fixture: REPEAT over codegen expressions (rates to Jacobian) with CSE
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

def rates():
    k = {}
    # $$JAFF REPEAT idx, rate IN rates
    k[0] = 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)
    k[1] = 2.6e-13*(tgas/10000.0)**(-0.8)
    k[2] = nden[3]
    # $$JAFF END
    return k


def rates_cse():
    k = {}
    # $$JAFF REPEAT idx, rate, cse IN rates
    k[0] = 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)
    k[1] = 4.1207223003989e-10*tgas**(-0.8)
    k[2] = nden[3]
    # $$JAFF END
    return k


def rates_inline():
    # $$JAFF REPEAT rate IN rates
    return [5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas), 2.6e-13*(tgas/10000.0)**(-0.8), nden[3]]
    # $$JAFF END


def fluxes():
    # $$JAFF REPEAT idx, flux_expression IN flux_expressions
    flux[0] = k[0] * y[idx_h] * y[idx_e]
    flux[1] = k[1] * y[idx_hj] * y[idx_e]
    flux[2] = k[2] * y[idx_hejj] * y[idx_e]
    # $$JAFF END


def ode_expressions():
    # $$JAFF REPEAT idx, ode_expression IN ode_expressions
    dy[0] =  - flux[0] + flux[1]
    dy[1] =  - flux[0] + flux[0] + flux[0] - flux[1] - flux[2]
    dy[2] =  + flux[0] - flux[1]
    dy[3] =  - flux[2]
    dy[4] =  + flux[2]
    # $$JAFF END


def odes():
    # $$JAFF REPEAT idx, ode IN odes
    dy[0] = 4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] - 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1]
    dy[1] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1] - nden[1]*nden[3]**2
    dy[2] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1]
    dy[3] = -nden[1]*nden[3]**2
    dy[4] = nden[1]*nden[3]**2
    # $$JAFF END


def odes_cse():
    # $$JAFF REPEAT idx, ode, cse IN odes
    tmp0_v = 4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] - 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1]
    tmp1_v = nden[1]*nden[3]**2
    dy[0] = tmp0_v
    dy[1] = -tmp0_v - tmp1_v
    dy[2] = -tmp0_v
    dy[3] = -tmp1_v
    dy[4] = tmp1_v
    # $$JAFF END


def rhs():
    # $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
    f[0] = 4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] - 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1]
    f[1] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1] - nden[1]*nden[3]**2
    f[2] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1]
    f[3] = -nden[1]*nden[3]**2
    f[4] = nden[1]*nden[3]**2
    f[5] = 0
    # $$JAFF END


def rhs_cse_replace():
    # $$JAFF REPEAT idx, rhs, cse IN rhses {
        # REPLACE {
            # "nden\[(\d+)\]" -> "y[\1]"
            # "tgas" -> "T"
        # }
        # THERMAL { MODE -> dedt }
    # }
    x0 = 4.1207223003989e-10*T**(-0.8)*y[1]*y[2] - 5.85e-11*math.sqrt(T)*math.exp(-157809.1/T)*y[0]*y[1]
    x1 = y[1]*y[3]**2
    f[1] = x0
    f[2] = -x0 - x1
    f[3] = -x0
    f[4] = -x1
    f[5] = x1
    f[6] = 0
    # $$JAFF END


def jacobian():
    # $$JAFF REPEAT idx, expr IN jacobian
    J[0][0] = -5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[1]
    J[0][1] = 4.1207223003989e-10*tgas**(-0.8)*nden[2] - 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]
    J[0][2] = 4.1207223003989e-10*tgas**(-0.8)*nden[1]
    J[1][0] = 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[1]
    J[1][1] = -4.1207223003989e-10*tgas**(-0.8)*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0] - nden[3]**2
    J[1][2] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]
    J[1][3] = -2*nden[1]*nden[3]
    J[2][0] = 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[1]
    J[2][1] = -4.1207223003989e-10*tgas**(-0.8)*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]
    J[2][2] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]
    J[3][1] = -nden[3]**2
    J[3][3] = -2*nden[1]*nden[3]
    J[4][1] = nden[3]**2
    J[4][3] = 2*nden[1]*nden[3]
    # $$JAFF END


def jacobian_cse():
    # $$JAFF REPEAT idx, expr, cse IN jacobian
    c0 = 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)
    c1 = c0*nden[1]
    c2 = 4.1207223003989e-10*tgas**(-0.8)
    c3 = -c0*nden[0] + c2*nden[2]
    c4 = c2*nden[1]
    c5 = nden[3]**2
    c6 = -c4
    c7 = 2*nden[1]*nden[3]
    c8 = -c7
    J[0, 0] = -c1
    J[0, 1] = c3
    J[0, 2] = c4
    J[1, 0] = c1
    J[1, 1] = -c3 - c5
    J[1, 2] = c6
    J[1, 3] = c8
    J[2, 0] = c1
    J[2, 1] = -c3
    J[2, 2] = c6
    J[3, 1] = -c5
    J[3, 3] = c8
    J[4, 1] = c5
    J[4, 3] = c7
    J1[1, 1] = -c1
    J1[1, 2] = c3
    J1[1, 3] = c4
    J1[2, 1] = c1
    J1[2, 2] = -c3 - c5
    J1[2, 3] = c6
    J1[2, 4] = c8
    J1[3, 1] = c1
    J1[3, 2] = -c3
    J1[3, 3] = c6
    J1[4, 2] = -c5
    J1[4, 4] = c8
    J1[5, 2] = c5
    J1[5, 4] = c7
    # $$JAFF END
