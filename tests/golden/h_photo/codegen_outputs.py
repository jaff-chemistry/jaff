# ABOUTME: jaffgen test template exposing every expression generator as a Python function
# ABOUTME: Rendered per network, then evaluated against tests/golden with shared inputs
import math

# $JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign $[POS j NEG k REPLACE idx_ek idx_e]$
idx_h = 0
idx_hj = 1
idx_e = 2
# $JAFF END

SPECIE_MASSES = {}
# $JAFF REPEAT idx, specie_mass IN specie_masses
SPECIE_MASSES[0] = 1.673773e-24
SPECIE_MASSES[1] = 1.6728620616289998e-24
SPECIE_MASSES[2] = 9.10938371e-28
# $JAFF END


def rates():
    k = {}
    # $JAFF REPEAT idx, rate, cse IN rates

    k[0] = 1.68692851322685e-18*c_hat*photden[0]
    k[1] = 1.65941781598291e-10*tgas**(-0.7)
    # $JAFF END
    return k


def fluxes():
    k = rates()
    y = nden
    out = {}
    # $JAFF REPEAT idx, flux_expression IN flux_expressions
    out[0] = k[0] * y[idx_h]
    out[1] = k[1] * y[idx_hj] * y[idx_e]
    # $JAFF END
    return out


def ode_expressions():
    flux = fluxes()
    out = {}
    # $JAFF REPEAT idx, ode_expression IN ode_expressions
    out[0] =  - flux[0] + flux[1]
    out[1] =  + flux[0] - flux[1]
    out[2] =  + flux[0] - flux[1]
    # $JAFF END
    return out


def odes():
    out = {}
    # $JAFF REPEAT idx, ode, cse IN odes
    cse0 = 1.68692851322685e-18*c_hat*nden[0]*photden[0] - 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2]

    out[0] = -cse0
    out[1] = cse0
    out[2] = cse0
    # $JAFF END
    return out


def radodes():
    out = {}
    # $JAFF REPEAT idx, radode, cse IN radodes
    cse0 = 1.68692851322685e-18*c_hat*nden[0]

    out[0] = -cse0*photden[0]
    out[1] = -cse0*rflux[0]
    # $JAFF END
    return out


def dedt():
    # $JAFF SUB dedt
    return {0: 1.07963424846518e-29*c_hat*nden[0]*photden[0] - 2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2]}
    # $JAFF END


def rhs_volumetric():
    out = {}
    # $JAFF REPEAT idx, rhs, cse IN rhses $[RADIATION True]$
    cse0 = c_hat*nden[0]*photden[0]
    cse1 = 1.68692851322685e-18*cse0
    cse2 = cse1 - 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2]

    out[0] = -cse2
    out[1] = cse2
    out[2] = cse2
    out[3] = 1.07963424846518e-29*cse0 - 2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2]
    out[4] = -cse1
    out[5] = -1.68692851322685e-18*c_hat*nden[0]*rflux[0]
    # $JAFF END
    return out


def rhs_specific_mass():
    out = {}
    # $JAFF REPEAT idx, rhs, cse IN rhses $[RADIATION True DEDT_TYPE specific]$
    cse0 = c_hat*nden[0]*photden[0]
    cse1 = 1.68692851322685e-18*cse0
    cse2 = cse1 - 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2]

    out[0] = -cse2
    out[1] = cse2
    out[2] = cse2
    out[3] = (1.07963424846518e-29*cse0 - 2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2])/(1.673773e-24*nden[0] + 1.672862061629e-24*nden[1] + 9.10938371e-28*nden[2])
    out[4] = -cse1
    out[5] = -1.68692851322685e-18*c_hat*nden[0]*rflux[0]
    # $JAFF END
    return out


def rhs_specific_number():
    out = {}
    # $JAFF REPEAT idx, rhs, cse IN rhses $[RADIATION True DEDT_TYPE per_particle]$
    cse0 = c_hat*nden[0]*photden[0]
    cse1 = 1.68692851322685e-18*cse0
    cse2 = cse1 - 1.65941781598291e-10*tgas**(-0.7)*nden[1]*nden[2]
    cse3 = nden[0] + nden[1] + nden[2]

    out[0] = -cse2
    out[1] = cse2
    out[2] = cse2
    out[3] = -1.49999999999992*cse2*tgas*(1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2])/cse3**2 + (1.07963424846518e-29*cse0 - 2.29107354821899e-26*tgas**0.3*(0.684 - 0.0416*math.log(0.0001*tgas))*nden[1]*nden[2])/cse3
    out[4] = -cse1
    out[5] = -1.68692851322685e-18*c_hat*nden[0]*rflux[0]
    # $JAFF END
    return out


def jacobian_species():
    out = {}
    # $JAFF REPEAT idx, expr, cse IN jacobian $[RADIATION True]$
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
    # $JAFF END
    return out


def jacobian_volumetric():
    out = {}
    # $JAFF REPEAT idx, expr, cse IN jacobian $[RADIATION True THERMAL dedt]$
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
    # $JAFF END
    return out


def jacobian_specific_mass():
    out = {}
    # $JAFF REPEAT idx, expr, cse IN jacobian $[RADIATION True THERMAL dedt DEDT_TYPE specific]$
    cse0 = 1.68692851322685e-18*c_hat
    cse1 = cse0*photden[0]
    cse2 = tgas**(-1.7)
    cse3 = 1.673773e-24*nden[0] + 1.672862061629e-24*nden[1] + 9.10938371e-28*nden[2]
    cse4 = 1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2]
    cse5 = 1/cse4
    cse6 = 2.0709734999999e-16*tgas/cse3
    cse7 = cse4*tgas/cse3**2
    cse8 = cse6 - 2.51065949999987e-24*cse7
    cse9 = cse1 - 7.74394980792062e-11*cse2*cse3*cse5*cse8*nden[1]*nden[2]
    cse10 = tgas**(-0.7)
    cse11 = 1.65941781598291e-10*cse10
    cse12 = cse6 - 2.50929309244337e-24*cse7
    cse13 = nden[1]*nden[2]
    cse14 = cse3*cse5
    cse15 = 7.74394980792062e-11*cse13*cse14*cse2
    cse16 = cse11*nden[2] + cse12*cse15
    cse17 = cse6 - 1.36640755649993e-27*cse7
    cse18 = cse11*nden[1] + cse15*cse17
    cse19 = cse0*nden[0]
    cse20 = -cse19
    cse21 = -cse16
    cse22 = -cse18
    cse23 = 1.673773e-24*nden[0] + 1.672862061629e-24*nden[1] + 9.10938371e-28*nden[2]
    cse24 = 1/cse23
    cse25 = 0.684 - 0.0416*math.log(0.0001*tgas)
    cse26 = 2.29107354821899e-26*cse25*tgas**0.3
    cse27 = (1.07963424846518e-29*c_hat*photden[0]*nden[0] - cse13*cse26)/cse23**2
    cse28 = 0.6666666666667*cse14*cse24*(-6.87322064465696e-27*cse10*cse13*cse25 + 9.53086596059098e-28*cse10*nden[1]*nden[2])
    cse29 = cse24*cse26

    out[(0, 0)] = -cse9
    out[(0, 1)] = cse16
    out[(0, 2)] = cse18
    out[(0, 3)] = -cse15
    out[(0, 4)] = cse20
    out[(1, 0)] = cse9
    out[(1, 1)] = cse21
    out[(1, 2)] = cse22
    out[(1, 3)] = cse15
    out[(1, 4)] = cse19
    out[(2, 0)] = cse9
    out[(2, 1)] = cse21
    out[(2, 2)] = cse22
    out[(2, 3)] = cse15
    out[(2, 4)] = cse19
    out[(3, 0)] = 1.07963424846518e-29*c_hat*cse24*photden[0] - 1.673773e-24*cse27 - cse28*cse8
    out[(3, 1)] = -cse12*cse28 - 1.672862061629e-24*cse27 - cse29*nden[2]
    out[(3, 2)] = -cse17*cse28 - 9.10938371e-28*cse27 - cse29*nden[1]
    out[(3, 3)] = cse28
    out[(3, 4)] = 1.07963424846518e-29*c_hat*cse24*nden[0]
    out[(4, 0)] = -cse1
    out[(4, 4)] = cse20
    out[(5, 0)] = -cse0*rflux[0]
    out[(5, 5)] = cse20
    # $JAFF END
    return out


def jacobian_specific_number():
    out = {}
    # $JAFF REPEAT idx, expr, cse IN jacobian $[RADIATION True THERMAL dedt DEDT_TYPE per_particle]$
    cse0 = 1.68692851322685e-18*c_hat
    cse1 = cse0*photden[0]
    cse2 = nden[0] + nden[1] + nden[2]
    cse3 = 1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2]
    cse4 = 2.0709734999999e-16*tgas/cse2 - 1.49999999999992*cse3*tgas/cse2**2
    cse5 = nden[1]*nden[2]
    cse6 = cse2/cse3
    cse7 = 7.74394980792062e-11*cse5*cse6*tgas**(-1.7)
    cse8 = cse4*cse7
    cse9 = cse1 - cse8
    cse10 = tgas**(-0.7)
    cse11 = 1.65941781598291e-10*cse10
    cse12 = cse11*nden[2]
    cse13 = cse12 + cse8
    cse14 = cse11*nden[1] + cse8
    cse15 = cse0*nden[0]
    cse16 = -cse15
    cse17 = -cse13
    cse18 = -cse14
    cse19 = nden[0] + nden[1] + nden[2]
    cse20 = 1/cse19
    cse21 = c_hat*photden[0]
    cse22 = cse19**(-2)
    cse23 = 1.380649e-16*nden[0] + 1.380649e-16*nden[1] + 1.380649e-16*nden[2]
    cse24 = 2.53039276984015e-18*cse22*cse23*tgas
    cse25 = tgas**0.3
    cse26 = 0.684 - 0.0416*math.log(0.0001*tgas)
    cse27 = 2.29107354821899e-26*cse25*cse26
    cse28 = cse1*nden[0] - cse12*nden[1]
    cse29 = cse22*cse28
    cse30 = cse10*cse5
    cse31 = 0.6666666666667*cse6*(cse20*(9.53086596059098e-28*cse10*nden[1]*nden[2] - 6.87322064465696e-27*cse26*cse30) - 1.74238870678197e-10*cse22*cse23*cse30 - 1.49999999999992*cse23*cse29)
    cse32 = cse22*(1.07963424846518e-29*cse21*nden[0] - cse27*cse5) + 2.0709734999999e-16*cse29*tgas + cse31*cse4 - 2.99999999999985*cse23*cse28*tgas/cse19**3
    cse33 = cse20*cse27

    out[(0, 0)] = -cse9
    out[(0, 1)] = cse13
    out[(0, 2)] = cse14
    out[(0, 3)] = -cse7
    out[(0, 4)] = cse16
    out[(1, 0)] = cse9
    out[(1, 1)] = cse17
    out[(1, 2)] = cse18
    out[(1, 3)] = cse7
    out[(1, 4)] = cse15
    out[(2, 0)] = cse9
    out[(2, 1)] = cse17
    out[(2, 2)] = cse18
    out[(2, 3)] = cse7
    out[(2, 4)] = cse15
    out[(3, 0)] = 1.07963424846518e-29*c_hat*cse20*photden[0] - cse21*cse24 - cse32
    out[(3, 1)] = 2.48912672397424e-10*cse22*cse23*cse25*nden[2] - cse32 - cse33*nden[2]
    out[(3, 2)] = 2.48912672397424e-10*cse22*cse23*cse25*nden[1] - cse32 - cse33*nden[1]
    out[(3, 3)] = cse31
    out[(3, 4)] = 1.07963424846518e-29*c_hat*cse20*nden[0] - c_hat*cse24*nden[0]
    out[(4, 0)] = -cse1
    out[(4, 4)] = cse16
    out[(5, 0)] = -cse0*rflux[0]
    out[(5, 5)] = cse16
    # $JAFF END
    return out
