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
    # $$JAFF REPEAT idx, rate, cse IN rates
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
    dy[0] = -nden[0]**2*nden[1]
    dy[1] = -nden[0]**2*nden[1]
    dy[2] = nden[0]**2*nden[1]
    # $$JAFF END


def odes_cse():
    # $$JAFF REPEAT idx, ode, cse IN odes
    tmp0_v = nden[0]**2*nden[1]
    tmp1_v = -tmp0_v
    dy[0] = tmp1_v
    dy[1] = tmp1_v
    dy[2] = tmp0_v
    # $$JAFF END


def rhs():
    # $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
    f[0] = -nden[0]**2*nden[1]
    f[1] = -nden[0]**2*nden[1]
    f[2] = nden[0]**2*nden[1]
    f[3] = nden[0]**3*nden[1]
    # $$JAFF END


def rhs_cse_replace():
    # $$JAFF REPEAT idx, rhs, cse IN rhses {
        # REPLACE {
            # "nden\[(\d+)\]" -> "y[\1]"
            # "tgas" -> "T"
        # }
        # THERMAL { MODE -> dedt }
    # }
    x0 = y[0]**2*y[1]
    x1 = -x0
    f[1] = x1
    f[2] = x1
    f[3] = x0
    f[4] = y[0]**3*y[1]
    # $$JAFF END


def jacobian():
    # $$JAFF REPEAT idx, expr IN jacobian
    J[0][0] = -2*nden[0]*nden[1]
    J[0][1] = -nden[0]**2
    J[1][0] = -2*nden[0]*nden[1]
    J[1][1] = -nden[0]**2
    J[2][0] = 2*nden[0]*nden[1]
    J[2][1] = nden[0]**2
    # $$JAFF END


def jacobian_cse():
    # $$JAFF REPEAT idx, expr, cse IN jacobian
    c0 = 2*nden[0]*nden[1]
    c1 = -c0
    c2 = nden[0]**2
    c3 = -c2
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
