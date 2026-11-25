# ABOUTME: Template-engine regression fixture: THERMAL options and the dedt/dtdt SUB tokens (energy network)
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

# $$JAFF SUB dedt, dtdt
dedt = nden[0]**3*nden[1]
dtdt = (2.0709734999999e-16*tgas*nden[0]**2*nden[1] + nden[0]**3*nden[1])/(2.0709734999999e-16*nden[0] + 2.0709734999999e-16*nden[1] + 2.0709734999999e-16*nden[2])
# $$JAFF END

# $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> none } }
f_none[0] = -nden[1]*nden[0]*nden[0]
f_none[1] = -nden[1]*nden[0]*nden[0]
f_none[2] = nden[1]*(nden[0]*nden[0])
# $$JAFF END

# $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
f_dedt[0] = -nden[1]*nden[0]*nden[0]
f_dedt[1] = -nden[1]*nden[0]*nden[0]
f_dedt[2] = nden[1]*(nden[0]*nden[0])
f_dedt[3] = nden[1]*((nden[0]*nden[0])*nden[0])
# $$JAFF END

# $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dtdt } }
f_dtdt[0] = -nden[1]*nden[0]*nden[0]
f_dtdt[1] = -nden[1]*nden[0]*nden[0]
f_dtdt[2] = nden[1]*(nden[0]*nden[0])
f_dtdt[3] = (2.0709734999999e-16*tgas*nden[1]*(nden[0]*nden[0]) + nden[1]*((nden[0]*nden[0])*nden[0]))/(2.0709734999999e-16*nden[0] + 2.0709734999999e-16*nden[1] + 2.0709734999999e-16*nden[2])
# $$JAFF END

# $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt  ENERGY -> specific } }
f_specific[0] = -nden[1]*nden[0]*nden[0]
f_specific[1] = -nden[1]*nden[0]*nden[0]
f_specific[2] = nden[1]*(nden[0]*nden[0])
f_specific[3] = -2.44367099999988e-24*tgas*(1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2])*nden[1]*nden[0]*nden[0]/(1.673773e-24*nden[0] + 3.343586e-24*nden[1] + 6.646473e-24*nden[2])**2 + nden[1]*((nden[0]*nden[0])*nden[0])/(1.673773e-24*nden[0] + 3.343586e-24*nden[1] + 6.646473e-24*nden[2])
# $$JAFF END

# $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt  ENERGY -> per_particle } }
f_per_particle[0] = -nden[1]*nden[0]*nden[0]
f_per_particle[1] = -nden[1]*nden[0]*nden[0]
f_per_particle[2] = nden[1]*(nden[0]*nden[0])
f_per_particle[3] = 1.49999999999992*tgas*(1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2])*nden[1]*(nden[0]*nden[0])/(nden[0] + nden[1] + nden[2])**2 + nden[1]*((nden[0]*nden[0])*nden[0])/(nden[0] + nden[1] + nden[2])
# $$JAFF END

# $$JAFF REPEAT idx, cse IN cses {
    # REPLACE { "nden\[(\d+)\]" -> "y[\1]" }
    # GROUP { NAME -> jacobian }
# }
c0 = y[0]*y[1]
c1 = 2*c0
c2 = -c1
c3 = y[0]*y[0]
c4 = -c3
c5 = 2.0709734999999e-16*y[1]
c6 = c5 + 2.0709734999999e-16*y[0] + 2.0709734999999e-16*y[2]
c7 = 1/c6
c8 = (y[0]*y[0])*y[0]
c9 = c3*tgas
c10 = -2.0709734999999e-16*(c5*c9 + c8*y[1])*1/(c6*c6)
# $$JAFF END
# $$JAFF REPEAT idx, jac IN jacobian {
    # REPLACE { "nden\[(\d+)\]" -> "y[\1]" }
    # THERMAL { MODE -> dtdt  ENERGY -> specific }
    # CSE { GROUP -> jacobian }
# }
J[0, 0] = c2
J[0, 1] = c4
J[1, 0] = c2
J[1, 1] = c4
J[2, 0] = c1
J[2, 1] = c3
J[3, 0] = c10 + c7*(4.14194699999979e-16*c0*tgas + 3*c3*y[1])
J[3, 1] = c10 + c7*(c8 + 2.0709734999999e-16*c9)
J[3, 2] = c10
J[3, 3] = c3*c5*c7
# $$JAFF END

# $$JAFF REPEAT idx, jac IN jacobian { THERMAL { MODE -> none } }
J_none[0, 0] = -2*nden[0]*nden[1]
J_none[0, 1] = -(nden[0]*nden[0])
J_none[1, 0] = -2*nden[0]*nden[1]
J_none[1, 1] = -(nden[0]*nden[0])
J_none[2, 0] = 2*nden[0]*nden[1]
J_none[2, 1] = nden[0]*nden[0]
# $$JAFF END
