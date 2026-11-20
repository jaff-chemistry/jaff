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
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    # $$JAFF END
    # $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
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
    dy[1] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1] - nden[1]*nden[3]*nden[3]
    dy[2] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1]
    dy[3] = -nden[1]*nden[3]*nden[3]
    dy[4] = nden[1]*(nden[3]*nden[3])
    # $$JAFF END


def odes_cse():
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> odes } }
    tmp0_v = 1/tgas
    tmp3_v = math.log(tgas)
    tmp4_v = math.exp(0.2*tmp3_v)
    tmp5_v = math.sqrt(tgas)
    tmp1_v = 4.1207223003989e-10*tmp0_v*tmp4_v*nden[1]*nden[2] - 5.85e-11*tmp5_v*math.exp(-157809.1*tmp0_v)*nden[0]*nden[1]
    tmp2_v = nden[1]*(nden[3]*nden[3])
    # $$JAFF END
    # $$JAFF REPEAT idx, ode IN odes { CSE { GROUP -> odes } }
    dy[0] = tmp1_v
    dy[1] = -tmp1_v - tmp2_v
    dy[2] = -tmp1_v
    dy[3] = -tmp2_v
    dy[4] = tmp2_v
    # $$JAFF END


def rhs():
    # $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
    f[0] = 4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] - 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1]
    f[1] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1] - nden[1]*nden[3]*nden[3]
    f[2] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]*nden[1]
    f[3] = -nden[1]*nden[3]*nden[3]
    f[4] = nden[1]*(nden[3]*nden[3])
    f[5] = 0
    # $$JAFF END


def rhs_cse_replace():
    # $$JAFF REPEAT idx, cse IN cses {
        # REPLACE {
            # "nden\[(\d+)\]" -> "y[\1]"
            # "tgas" -> "T"
        # }
        # GROUP { NAME -> rhs }
    # }
    x0 = 1/T
    x3 = math.log(T)
    x4 = math.exp(0.2*x3)
    x5 = math.sqrt(T)
    x1 = 4.1207223003989e-10*x0*x4*y[1]*y[2] - 5.85e-11*x5*math.exp(-157809.1*x0)*y[0]*y[1]
    x2 = y[1]*(y[3]*y[3])
    # $$JAFF END
    # $$JAFF REPEAT idx, rhs IN rhses {
        # REPLACE {
            # "nden\[(\d+)\]" -> "y[\1]"
            # "tgas" -> "T"
        # }
        # THERMAL { MODE -> dedt }
        # CSE { GROUP -> rhs }
    # }
    f[1] = x1
    f[2] = -x1 - x2
    f[3] = -x1
    f[4] = -x2
    f[5] = x2
    f[6] = 0
    # $$JAFF END


def jacobian():
    # $$JAFF REPEAT idx, jac IN jacobian
    J[0][0] = -5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[1]
    J[0][1] = 4.1207223003989e-10*tgas**(-0.8)*nden[2] - 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]
    J[0][2] = 4.1207223003989e-10*tgas**(-0.8)*nden[1]
    J[1][0] = 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[1]
    J[1][1] = -4.1207223003989e-10*tgas**(-0.8)*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0] - (nden[3]*nden[3])
    J[1][2] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]
    J[1][3] = -2*nden[1]*nden[3]
    J[2][0] = 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[1]
    J[2][1] = -4.1207223003989e-10*tgas**(-0.8)*nden[2] + 5.85e-11*math.sqrt(tgas)*math.exp(-157809.1/tgas)*nden[0]
    J[2][2] = -4.1207223003989e-10*tgas**(-0.8)*nden[1]
    J[3][1] = -(nden[3]*nden[3])
    J[3][3] = -2*nden[1]*nden[3]
    J[4][1] = nden[3]*nden[3]
    J[4][3] = 2*nden[1]*nden[3]
    # $$JAFF END


def jacobian_cse():
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    c0 = 1/tgas
    c12 = math.sqrt(tgas)
    c1 = 5.85e-11*c12*math.exp(-157809.1*c0)
    c2 = c1*nden[1]
    c13 = math.log(tgas)
    c14 = math.exp(0.2*c13)
    c3 = 4.1207223003989e-10*c0*c14
    c4 = c3*nden[2]
    c5 = c1*nden[0]
    c6 = -c4 + c5
    c7 = c3*nden[1]
    c8 = nden[3]*nden[3]
    c9 = -c7
    c10 = 2*nden[1]*nden[3]
    c11 = -c10
    # $$JAFF END
    # $$JAFF REPEAT idx, jac IN jacobian { CSE { GROUP -> jacobian } }
    J[0, 0] = -c2
    J[0, 1] = -c6
    J[0, 2] = c7
    J[1, 0] = c2
    J[1, 1] = -c4 + c5 - c8
    J[1, 2] = c9
    J[1, 3] = c11
    J[2, 0] = c2
    J[2, 1] = c6
    J[2, 2] = c9
    J[3, 1] = -c8
    J[3, 3] = c11
    J[4, 1] = c8
    J[4, 3] = c10
    J1[1, 1] = -c2
    J1[1, 2] = -c6
    J1[1, 3] = c7
    J1[2, 1] = c2
    J1[2, 2] = -c4 + c5 - c8
    J1[2, 3] = c9
    J1[2, 4] = c11
    J1[3, 1] = c2
    J1[3, 2] = c6
    J1[3, 3] = c9
    J1[4, 2] = -c8
    J1[4, 4] = c11
    J1[5, 2] = c8
    J1[5, 4] = c10
    # $$JAFF END
