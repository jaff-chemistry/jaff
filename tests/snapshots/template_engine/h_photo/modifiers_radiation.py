# ABOUTME: Template-engine regression fixture: radodes and RADIATION options (radiation network)
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

# $$JAFF REPEAT idx, radode IN radodes
drad[0] = -5.05728445450563e-8*nden[0]*photden[0]
drad[1] = -5.05728445450563e-8*nden[0]*rflux[0]
# $$JAFF END

# $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> radodes } }
r0 = 5.05728445450563e-8*nden[0]
# $$JAFF END
# $$JAFF REPEAT idx, radode IN radodes { CSE { GROUP -> radodes } }
drad[0] = -r0*photden[0]
drad[1] = -r0*rflux[0]
# $$JAFF END

# $$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> True }  THERMAL { MODE -> dedt } }
f_rad[0] = 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] - 5.05728445450563e-8*nden[0]*photden[0]
f_rad[1] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
f_rad[2] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
f_rad[3] = -2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2] + 3.2366620508836e-19*nden[0]*photden[0]
f_rad[4] = -5.05728445450563e-8*nden[0]*photden[0]
f_rad[5] = -5.05728445450563e-8*nden[0]*rflux[0]
# $$JAFF END

# $$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> FALSE }  THERMAL { MODE -> dedt } }
f_norad[0] = 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] - 5.05728445450563e-8*nden[0]*photden[0]
f_norad[1] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
f_norad[2] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
f_norad[3] = -2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2] + 3.2366620508836e-19*nden[0]*photden[0]
# $$JAFF END

# $$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> true  ORDER -> 1 }  THERMAL { MODE -> dedt } }
f_rad1[0] = 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] - 5.05728445450563e-8*nden[0]*photden[0]
f_rad1[1] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
f_rad1[2] = -1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2] + 5.05728445450563e-8*nden[0]*photden[0]
f_rad1[3] = -2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2] + 3.2366620508836e-19*nden[0]*photden[0]
f_rad1[4] = -5.05728445450563e-8*nden[0]*rflux[0]
f_rad1[5] = -5.05728445450563e-8*nden[0]*photden[0]
# $$JAFF END

# $$JAFF REPEAT idx, cse IN cses {
    # REPLACE {
        # "photden\[\s*(\d+)\s*\]" -> "rn[2*\1]"
        # "rflux\[\s*(\d+)\s*\]" -> "rn[2*\1+1]"
        # "nden\[\s*(\d+)\s*\]" -> "xn[\1]"
        # "tgas" -> "T"
        # "c_hat" -> "chat"
    # }
    # GROUP { NAME -> jacobian }
# }
c0 = 5.05728445450563e-8*rn[2*0]
c1 = -c0
c2 = 1.65941781598291e-10*T**(-0.7)
c3 = c2*xn[2]
c4 = c2*xn[1]
c5 = 5.05728445450563e-8*xn[0]
c6 = -c5
c7 = -c3
c8 = -c4
# $$JAFF END
# $$JAFF REPEAT idx, jac IN jacobian {
    # REPLACE {
        # "photden\[\s*(\d+)\s*\]" -> "rn[2*\1]"
        # "rflux\[\s*(\d+)\s*\]" -> "rn[2*\1+1]"
        # "nden\[\s*(\d+)\s*\]" -> "xn[\1]"
        # "tgas" -> "T"
        # "c_hat" -> "chat"
    # }
    # RADIATION { ENABLED -> True }
    # CSE { GROUP -> jacobian }
# }
J[0, 0] = c1
J[0, 1] = c3
J[0, 2] = c4
J[0, 3] = c6
J[1, 0] = c0
J[1, 1] = c7
J[1, 2] = c8
J[1, 3] = c5
J[2, 0] = c0
J[2, 1] = c7
J[2, 2] = c8
J[2, 3] = c5
J[3, 0] = c1
J[3, 3] = c6
J[4, 0] = -5.05728445450563e-8*rn[2*0+1]
J[4, 4] = c6
# $$JAFF END
