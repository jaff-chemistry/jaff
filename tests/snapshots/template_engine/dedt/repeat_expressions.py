# ABOUTME: Template-engine regression fixture: REPEAT over codegen expressions (rates to Jacobian) with CSE
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

def rates():
    k = {}
    # $$JAFF REPEAT idx, rate IN rates
    k[0] = nden[0]
    # $$JAFF END
    return k


def rates_cse():
    k = {}
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    # $$JAFF END
    # $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
    k[0] = nden[0]
    # $$JAFF END
    return k


def rates_inline():
    # $$JAFF REPEAT rate IN rates
    return [nden[0]]
    # $$JAFF END


def fluxes():
    # $$JAFF REPEAT idx, flux_expression IN flux_expressions
    flux[0] = k[0] * y[idx_h] * y[idx_d]
    # $$JAFF END


def ode_expressions():
    # $$JAFF REPEAT idx, ode_expression IN ode_expressions
    dy[0] =  - flux[0]
    dy[1] =  - flux[0]
    dy[2] =  + flux[0]
    # $$JAFF END


def odes():
    # $$JAFF REPEAT idx, ode IN odes
    dy[0] = -nden[1]*nden[0]*nden[0]
    dy[1] = -nden[1]*nden[0]*nden[0]
    dy[2] = nden[1]*(nden[0]*nden[0])
    # $$JAFF END


def odes_cse():
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> odes } }
    tmp0_v = nden[1]*(nden[0]*nden[0])
    tmp1_v = -tmp0_v
    # $$JAFF END
    # $$JAFF REPEAT idx, ode IN odes { CSE { GROUP -> odes } }
    dy[0] = tmp1_v
    dy[1] = tmp1_v
    dy[2] = tmp0_v
    # $$JAFF END


def rhs():
    # $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
    f[0] = -nden[1]*nden[0]*nden[0]
    f[1] = -nden[1]*nden[0]*nden[0]
    f[2] = nden[1]*(nden[0]*nden[0])
    f[3] = nden[1]*((nden[0]*nden[0])*nden[0])
    # $$JAFF END


def rhs_cse_replace():
    # $$JAFF REPEAT idx, cse IN cses {
        # REPLACE {
            # "nden\[(\d+)\]" -> "y[\1]"
            # "tgas" -> "T"
        # }
        # GROUP { NAME -> rhs }
    # }
    x0 = y[1]*(y[0]*y[0])
    x1 = -x0
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
    f[2] = x1
    f[3] = x0
    f[4] = y[1]*((y[0]*y[0])*y[0])
    # $$JAFF END


def jacobian():
    # $$JAFF REPEAT idx, jac IN jacobian
    J[0][0] = -2*nden[0]*nden[1]
    J[0][1] = -(nden[0]*nden[0])
    J[1][0] = -2*nden[0]*nden[1]
    J[1][1] = -(nden[0]*nden[0])
    J[2][0] = 2*nden[0]*nden[1]
    J[2][1] = nden[0]*nden[0]
    # $$JAFF END


def jacobian_cse():
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    c0 = 2*nden[0]*nden[1]
    c1 = -c0
    c2 = nden[0]*nden[0]
    c3 = -c2
    # $$JAFF END
    # $$JAFF REPEAT idx, jac IN jacobian { CSE { GROUP -> jacobian } }
    J[0, 0] = c1
    J[0, 1] = c3
    J[1, 0] = c1
    J[1, 1] = c3
    J[2, 0] = c0
    J[2, 1] = c2
    J1[1, 1] = c1
    J1[1, 2] = c3
    J1[2, 1] = c1
    J1[2, 2] = c3
    J1[3, 1] = c0
    J1[3, 2] = c2
    # $$JAFF END
