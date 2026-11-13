# ABOUTME: Template-engine regression fixture: REPEAT over codegen expressions (rates to Jacobian) with CSE
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

def rates():
    k = {}
    # $$JAFF REPEAT idx, rate IN rates
    k[0] = 5.05728445450563e-8*photden[0]
    k[1] = 2.63e-13*(tgas/10000.0)**(-0.7)
    # $$JAFF END
    return k


def rates_cse():
    k = {}
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    # $$JAFF END
    # $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
    k[0] = 5.05728445450563e-8*photden[0]
    k[1] = 1.65941781598291e-10*tgas**(-0.7)
    # $$JAFF END
    return k


def rates_inline():
    # $$JAFF REPEAT rate IN rates
    return [5.05728445450563e-8*photden[0], 2.63e-13*(tgas/10000.0)**(-0.7)]
    # $$JAFF END


def fluxes():
    # $$JAFF REPEAT idx, flux_expression IN flux_expressions
    flux[0] = k[0] * y[idx_h]
    flux[1] = k[1] * y[idx_hj] * y[idx_e]
    # $$JAFF END


def ode_expressions():
    # $$JAFF REPEAT idx, ode_expression IN ode_expressions
    dy[0] =  - flux[0] + flux[1]
    dy[1] =  + flux[0] - flux[1]
    dy[2] =  + flux[0] - flux[1]
    # $$JAFF END


def odes():
    # $$JAFF REPEAT idx, ode IN odes
    dy[0] = 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] - 5.05728445450563e-8*nden[0]*photden[0]
    dy[1] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
    dy[2] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
    # $$JAFF END


def odes_cse():
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> odes } }
    tmp0_v = 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] - 5.05728445450563e-8*nden[0]*photden[0]
    tmp1_v = -tmp0_v
    # $$JAFF END
    # $$JAFF REPEAT idx, ode IN odes { CSE { GROUP -> odes } }
    dy[0] = tmp0_v
    dy[1] = tmp1_v
    dy[2] = tmp1_v
    # $$JAFF END


def rhs():
    # $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
    f[0] = 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] - 5.05728445450563e-8*nden[0]*photden[0]
    f[1] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
    f[2] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
    f[3] = -2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2] + 3.2366620508836e-19*nden[0]*photden[0]
    # $$JAFF END


def rhs_cse_replace():
    # $$JAFF REPEAT idx, cse IN cses {
        # REPLACE {
            # "nden\[(\d+)\]" -> "y[\1]"
            # "tgas" -> "T"
        # }
        # GROUP { NAME -> rhs }
    # }
    x0 = y[1]*y[2]
    x1 = 1.65941781598291e-10*T**(-0.7)*x0 - 5.05728445450563e-8*y[0]*photden[0]
    x2 = -x1
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
    f[2] = x2
    f[3] = x2
    f[4] = -2.29107354821899e-26*T**0.3*x0*(0.684 - 0.0416*math.log(0.0001*T)) + 3.2366620508836e-19*y[0]*photden[0]
    # $$JAFF END


def jacobian():
    # $$JAFF REPEAT idx, jac IN jacobian
    J[0][0] = -5.05728445450563e-8*photden[0]
    J[0][1] = 1.65941781598291e-10*tgas**(-0.7)*nden[2]
    J[0][2] = 1.65941781598291e-10*tgas**(-0.7)*nden[1]
    J[1][0] = 5.05728445450563e-8*photden[0]
    J[1][1] = -1.65941781598291e-10*tgas**(-0.7)*nden[2]
    J[1][2] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]
    J[2][0] = 5.05728445450563e-8*photden[0]
    J[2][1] = -1.65941781598291e-10*tgas**(-0.7)*nden[2]
    J[2][2] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]
    # $$JAFF END


def jacobian_cse():
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    c0 = 5.05728445450563e-8*photden[0]
    c1 = 1.65941781598291e-10*tgas**(-0.7)
    c2 = c1*nden[2]
    c3 = c1*nden[1]
    c4 = -c2
    c5 = -c3
    # $$JAFF END
    # $$JAFF REPEAT idx, jac IN jacobian { CSE { GROUP -> jacobian } }
    J[0, 0] = -c0
    J[0, 1] = c2
    J[0, 2] = c3
    J[1, 0] = c0
    J[1, 1] = c4
    J[1, 2] = c5
    J[2, 0] = c0
    J[2, 1] = c4
    J[2, 2] = c5
    J1[1, 1] = -c0
    J1[1, 2] = c2
    J1[1, 3] = c3
    J1[2, 1] = c0
    J1[2, 2] = c4
    J1[2, 3] = c5
    J1[3, 1] = c0
    J1[3, 2] = c4
    J1[3, 3] = c5
    # $$JAFF END
