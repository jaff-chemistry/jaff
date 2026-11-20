# ABOUTME: jaffgen test template exposing every expression generator as a Python function
# ABOUTME: Rendered per network, then evaluated against tests/golden with shared inputs
import math

# $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
    # REPLACE { "idx_ek" -> "idx_e" }
    # SIGN { POS -> j  NEG -> k }
# }
idx_h = 0
idx_hj = 1
idx_e = 2
# $$JAFF END

SPECIE_MASSES = {}
# $$JAFF REPEAT idx, specie_mass IN specie_masses
SPECIE_MASSES[0] = 1.673773e-24
SPECIE_MASSES[1] = 1.6728620616289998e-24
SPECIE_MASSES[2] = 9.10938371e-28
# $$JAFF END


def rates():
    k = {}
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    # $$JAFF END
    # $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
    k[0] = 1.68692851322685e-18*c_hat*photden[0]
    k[1] = 1.65941781598291e-10*tgas**(-0.7)
    # $$JAFF END
    return k


def fluxes():
    k = rates()
    y = nden
    out = {}
    # $$JAFF REPEAT idx, flux_expression IN flux_expressions
    out[0] = k[0] * y[idx_h]
    out[1] = k[1] * y[idx_hj] * y[idx_e]
    # $$JAFF END
    return out


def ode_expressions():
    flux = fluxes()
    out = {}
    # $$JAFF REPEAT idx, ode_expression IN ode_expressions
    out[0] =  - flux[0] + flux[1]
    out[1] =  + flux[0] - flux[1]
    out[2] =  + flux[0] - flux[1]
    # $$JAFF END
    return out


def odes():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> odes } }
    cse0 = 1.68692851322685e-18*c_hat*nden[0]*photden[0] - 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2]
    # $$JAFF END
    # $$JAFF REPEAT idx, ode IN odes { CSE { GROUP -> odes } }
    out[0] = -cse0
    out[1] = cse0
    out[2] = cse0
    # $$JAFF END
    return out


def radodes():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> radodes } }
    cse0 = 1.68692851322685e-18*c_hat*nden[0]
    # $$JAFF END
    # $$JAFF REPEAT idx, radode IN radodes { CSE { GROUP -> radodes } }
    out[0] = -cse0*photden[0]
    out[1] = -cse0*rflux[0]
    # $$JAFF END
    return out


def dedt():
    # $$JAFF SUB dedt
    return {0: 1.07963424846518e-29*c_hat*nden[0]*photden[0] - 2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2]}
    # $$JAFF END


def rhs_volumetric():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rhs } }
    cse0 = c_hat*nden[0]*photden[0]
    cse1 = 1.68692851322685e-18*cse0
    cse2 = cse1 - 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2]
    # $$JAFF END
    # $$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> True }  THERMAL { MODE -> dedt }  CSE { GROUP -> rhs } }
    out[0] = -cse2
    out[1] = cse2
    out[2] = cse2
    out[3] = 1.07963424846518e-29*cse0 - 2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2]
    out[4] = -cse1
    out[5] = -1.68692851322685e-18*c_hat*nden[0]*rflux[0]
    # $$JAFF END
    return out


def rhs_specific_mass():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses {
        # GROUP { NAME -> rhs_1 }
    # }
    cse0 = c_hat*nden[0]*photden[0]
    cse1 = 1.68692851322685e-18*cse0
    cse2 = cse1 - 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2]
    # $$JAFF END
    # $$JAFF REPEAT idx, rhs IN rhses {
        # RADIATION { ENABLED -> True }
        # THERMAL { MODE -> dedt  ENERGY -> specific }
        # CSE { GROUP -> rhs_1 }
    # }
    out[0] = -cse2
    out[1] = cse2
    out[2] = cse2
    out[3] = (1.07963424846518e-29*cse0 - 2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2])/(1.673773e-24*nden[0] + 1.672862061629e-24*nden[1] + 9.10938371e-28*nden[2])
    out[4] = -cse1
    out[5] = -1.68692851322685e-18*c_hat*nden[0]*rflux[0]
    # $$JAFF END
    return out


def rhs_specific_number():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses {
        # GROUP { NAME -> rhs_2 }
    # }
    cse0 = c_hat*nden[0]*photden[0]
    cse1 = 1.68692851322685e-18*cse0
    cse2 = cse1 - 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2]
    cse3 = nden[0] + nden[1] + nden[2]
    # $$JAFF END
    # $$JAFF REPEAT idx, rhs IN rhses {
        # RADIATION { ENABLED -> True }
        # THERMAL { MODE -> dedt  ENERGY -> per_particle }
        # CSE { GROUP -> rhs_2 }
    # }
    out[0] = -cse2
    out[1] = cse2
    out[2] = cse2
    out[3] = -1.49999999999992*cse2*tgas*(1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2])*1/(cse3*cse3) + (1.07963424846518e-29*cse0 - 2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2])/cse3
    out[4] = -cse1
    out[5] = -1.68692851322685e-18*c_hat*nden[0]*rflux[0]
    # $$JAFF END
    return out


def jacobian_species():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    cse0 = 1.68692851322685e-18*c_hat
    cse1 = cse0*photden[0]
    cse2 = -cse1
    cse3 = 1.65941781598291e-10*tgas**(-0.7)
    cse4 = cse3*nden[2]
    cse5 = cse3*nden[1]
    cse6 = cse0*nden[0]
    cse7 = -cse6
    cse8 = -cse4
    cse9 = -cse5
    # $$JAFF END
    # $$JAFF REPEAT idx, jac IN jacobian { RADIATION { ENABLED -> True }  CSE { GROUP -> jacobian } }
    out[(0, 0)] = cse2
    out[(0, 1)] = cse4
    out[(0, 2)] = cse5
    out[(0, 3)] = cse7
    out[(1, 0)] = cse1
    out[(1, 1)] = cse8
    out[(1, 2)] = cse9
    out[(1, 3)] = cse6
    out[(2, 0)] = cse1
    out[(2, 1)] = cse8
    out[(2, 2)] = cse9
    out[(2, 3)] = cse6
    out[(3, 0)] = cse2
    out[(3, 3)] = cse7
    out[(4, 0)] = -cse0*rflux[0]
    out[(4, 4)] = cse7
    # $$JAFF END
    return out


def jacobian_volumetric():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses {
        # GROUP { NAME -> jacobian_1 }
    # }
    cse0 = 1.68692851322685e-18*c_hat
    cse1 = cse0*photden[0]
    cse2 = 1/(2.0709734999999e-16*nden[0] + 2.0709734999999e-16*nden[1] + 2.0709734999999e-16*nden[2])
    cse3 = tgas**(-0.7)
    cse4 = nden[1]*nden[2]
    cse5 = cse3*cse4
    cse6 = 2.40562722562981e-26*cse2*cse5
    cse7 = cse1 - cse6
    cse8 = 1.65941781598291e-10*cse3
    cse9 = cse6 + cse8*nden[2]
    cse10 = cse6 + cse8*nden[1]
    cse11 = 1.16159247118804e-10*cse2*cse4*tgas**(-1.7)
    cse12 = cse0*nden[0]
    cse13 = -cse12
    cse14 = -cse9
    cse15 = -cse10
    cse16 = 1.07963424846518e-29*c_hat
    cse17 = 0.684 - 0.0416*math.log(0.0001*tgas)
    cse18 = cse2*(-6.87322064465696e-27*cse17*cse5 + 9.53086596059098e-28*cse3*nden[1]*nden[2])
    cse19 = 2.0709734999999e-16*cse18*tgas
    cse20 = 2.29107354821899e-26*cse17*tgas**0.3
    # $$JAFF END
    # $$JAFF REPEAT idx, jac IN jacobian {
        # RADIATION { ENABLED -> True }
        # THERMAL { MODE -> dedt }
        # CSE { GROUP -> jacobian_1 }
    # }
    out[(0, 0)] = -cse7
    out[(0, 1)] = cse9
    out[(0, 2)] = cse10
    out[(0, 3)] = -cse11
    out[(0, 4)] = cse13
    out[(1, 0)] = cse7
    out[(1, 1)] = cse14
    out[(1, 2)] = cse15
    out[(1, 3)] = cse11
    out[(1, 4)] = cse12
    out[(2, 0)] = cse7
    out[(2, 1)] = cse14
    out[(2, 2)] = cse15
    out[(2, 3)] = cse11
    out[(2, 4)] = cse12
    out[(3, 0)] = cse16*photden[0] - cse19
    out[(3, 1)] = -cse19 - cse20*nden[2]
    out[(3, 2)] = -cse19 - cse20*nden[1]
    out[(3, 3)] = cse18
    out[(3, 4)] = cse16*nden[0]
    out[(4, 0)] = -cse1
    out[(4, 4)] = cse13
    out[(5, 0)] = -cse0*rflux[0]
    out[(5, 5)] = cse13
    # $$JAFF END
    return out


def jacobian_specific_mass():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses {
        # GROUP { NAME -> jacobian_2 }
    # }
    cse0 = 1.68692851322685e-18*c_hat
    cse1 = cse0*photden[0]
    cse2 = 1.673773e-24*nden[0] + 1.672862061629e-24*nden[1] + 9.10938371e-28*nden[2]
    cse3 = 1/cse2
    cse4 = 2.0709734999999e-16*cse3*tgas
    cse5 = 1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2]
    cse6 = 1/(cse2*cse2)
    cse7 = -cse4 + 2.51065949999987e-24*cse5*cse6*tgas
    cse8 = 1/cse5
    cse9 = nden[1]*nden[2]
    cse10 = 7.74394980792062e-11*cse2*cse8*cse9*tgas**(-1.7)
    cse11 = cse1 + cse10*cse7
    cse12 = tgas**(-0.7)
    cse13 = -cse4 + 2.50929309244337e-24*cse5*cse6*tgas
    cse14 = cse10*cse13 - 1.65941781598291e-10*cse12*nden[2]
    cse15 = -cse4 + 1.36640755649993e-27*cse5*cse6*tgas
    cse16 = cse10*cse15 - 1.65941781598291e-10*cse12*nden[1]
    cse17 = cse0*nden[0]
    cse18 = -cse17
    cse19 = 1.07963424846518e-29*c_hat*photden[0]
    cse20 = 0.684 - 0.0416*math.log(0.0001*tgas)
    cse21 = 2.29107354821899e-26*cse20*tgas**0.3
    cse22 = cse6*(cse19*nden[0] - cse21*cse9)
    cse23 = -6.87322064465696e-27*cse12*cse20*cse9 + 9.53086596059098e-28*cse12*nden[1]*nden[2]
    cse24 = 0.6666666666667*cse23*cse8
    cse25 = cse21*cse3
    # $$JAFF END
    # $$JAFF REPEAT idx, jac IN jacobian {
        # RADIATION { ENABLED -> True }
        # THERMAL { MODE -> dedt  ENERGY -> specific }
        # CSE { GROUP -> jacobian_2 }
    # }
    out[(0, 0)] = -cse11
    out[(0, 1)] = -cse14
    out[(0, 2)] = -cse16
    out[(0, 3)] = -cse10
    out[(0, 4)] = cse18
    out[(1, 0)] = cse11
    out[(1, 1)] = cse14
    out[(1, 2)] = cse16
    out[(1, 3)] = cse10
    out[(1, 4)] = cse17
    out[(2, 0)] = cse11
    out[(2, 1)] = cse14
    out[(2, 2)] = cse16
    out[(2, 3)] = cse10
    out[(2, 4)] = cse17
    out[(3, 0)] = cse19*cse3 - 1.673773e-24*cse22 + cse24*cse7
    out[(3, 1)] = 0.6666666666667*cse13*cse23*cse8 - 1.672862061629e-24*cse22 - cse25*nden[2]
    out[(3, 2)] = 0.6666666666667*cse15*cse23*cse8 - 9.10938371e-28*cse22 - cse25*nden[1]
    out[(3, 3)] = cse24
    out[(3, 4)] = 1.07963424846518e-29*c_hat*cse3*nden[0]
    out[(4, 0)] = -cse1
    out[(4, 4)] = cse18
    out[(5, 0)] = -cse0*rflux[0]
    out[(5, 5)] = cse18
    # $$JAFF END
    return out


def jacobian_specific_number():
    out = {}
    # $$JAFF REPEAT idx, cse IN cses {
        # GROUP { NAME -> jacobian_3 }
    # }
    cse0 = 1.68692851322685e-18*c_hat
    cse1 = cse0*photden[0]
    cse2 = nden[0] + nden[1] + nden[2]
    cse3 = 1/cse2
    cse4 = 2.0709734999999e-16*tgas
    cse5 = 1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2]
    cse6 = cse5*tgas
    cse7 = 1/(cse2*cse2)
    cse8 = 1.49999999999992*cse7
    cse9 = -cse3*cse4 + cse6*cse8
    cse10 = nden[1]*nden[2]
    cse11 = 1/cse5
    cse12 = cse11*cse2
    cse13 = 7.74394980792062e-11*cse10*cse12*tgas**(-1.7)
    cse14 = cse13*cse9
    cse15 = cse1 + cse14
    cse16 = tgas**(-0.7)
    cse17 = cse14 - 1.65941781598291e-10*cse16*nden[2]
    cse18 = cse14 - 1.65941781598291e-10*cse16*nden[1]
    cse19 = cse0*nden[0]
    cse20 = -cse19
    cse21 = c_hat*photden[0]
    cse22 = 2.53039276984015e-18*cse6*cse7
    cse23 = tgas**0.3
    cse24 = 0.684 - 0.0416*math.log(0.0001*tgas)
    cse25 = 2.29107354821899e-26*cse23*cse24
    cse26 = cse1*nden[0] - 1.65941781598291e-10*cse16*nden[1]*nden[2]
    cse27 = cse10*cse16
    cse28 = -cse26*cse5*cse8 - 1.74238870678197e-10*cse27*cse5*cse7 + cse3*(9.53086596059098e-28*cse16*nden[1]*nden[2] - 6.87322064465696e-27*cse24*cse27)
    cse29 = -0.6666666666667*cse11*cse2*cse28*cse9 + cse26*cse4*cse7 - 2.99999999999985*cse26*cse5*tgas*((cse2*cse2)*cse2)**(-1) + cse7*(-cse10*cse25 + 1.07963424846518e-29*cse21*nden[0])
    cse30 = cse25*cse3
    # $$JAFF END
    # $$JAFF REPEAT idx, jac IN jacobian {
        # RADIATION { ENABLED -> True }
        # THERMAL { MODE -> dedt  ENERGY -> per_particle }
        # CSE { GROUP -> jacobian_3 }
    # }
    out[(0, 0)] = -cse15
    out[(0, 1)] = -cse17
    out[(0, 2)] = -cse18
    out[(0, 3)] = -cse13
    out[(0, 4)] = cse20
    out[(1, 0)] = cse15
    out[(1, 1)] = cse17
    out[(1, 2)] = cse18
    out[(1, 3)] = cse13
    out[(1, 4)] = cse19
    out[(2, 0)] = cse15
    out[(2, 1)] = cse17
    out[(2, 2)] = cse18
    out[(2, 3)] = cse13
    out[(2, 4)] = cse19
    out[(3, 0)] = 1.07963424846518e-29*c_hat*cse3*photden[0] - cse21*cse22 - cse29
    out[(3, 1)] = 2.48912672397424e-10*cse23*cse5*cse7*nden[2] - cse29 - cse30*nden[2]
    out[(3, 2)] = 2.48912672397424e-10*cse23*cse5*cse7*nden[1] - cse29 - cse30*nden[1]
    out[(3, 3)] = 0.6666666666667*cse12*cse28
    out[(3, 4)] = -c_hat*cse22*nden[0] + 1.07963424846518e-29*c_hat*cse3*nden[0]
    out[(4, 0)] = -cse1
    out[(4, 4)] = cse20
    out[(5, 0)] = -cse0*rflux[0]
    out[(5, 5)] = cse20
    # $$JAFF END
    return out
