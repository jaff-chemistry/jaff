# ABOUTME: jaffgen test template: every auxiliary function of the GOW network
# ABOUTME: Rendered via GET aux_func; compared numerically against tests/golden
import math


def aux_functions():
    out = {}
    $$JAFF GET aux_func IN kcr_h_fac
    out["kcr_h_fac"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN qcr_h
    out["qcr_h"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN qcr_h2
    out["qcr_h2"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN qcr
    out["qcr"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN ncrh2
    out["ncrh2"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN kgr_gong
    out["kgr_gong"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN n2ncr
    out["n2ncr"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN fshield_h2
    out["fshield_h2"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN fshield_c
    out["fshield_c"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff0
    out["chemratecoeff0"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff1
    out["chemratecoeff1"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff4
    out["chemratecoeff4"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff8
    out["chemratecoeff8"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff10
    out["chemratecoeff10"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff13
    out["chemratecoeff13"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff14
    out["chemratecoeff14"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff15
    out["chemratecoeff15"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff16
    out["chemratecoeff16"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff17
    out["chemratecoeff17"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff18
    out["chemratecoeff18"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff19
    out["chemratecoeff19"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff20
    out["chemratecoeff20"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff21
    out["chemratecoeff21"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff22
    out["chemratecoeff22"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff23
    out["chemratecoeff23"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff27
    out["chemratecoeff27"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff28
    out["chemratecoeff28"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff32
    out["chemratecoeff32"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff35
    out["chemratecoeff35"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff37
    out["chemratecoeff37"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff39
    out["chemratecoeff39"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff40
    out["chemratecoeff40"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff41
    out["chemratecoeff41"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff42
    out["chemratecoeff42"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff48
    out["chemratecoeff48"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN chemratecoeff49
    out["chemratecoeff49"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN deltae0
    out["deltae0"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN deltae1
    out["deltae1"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN deltae2
    out["deltae2"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN deltae13
    out["deltae13"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN deltae14
    out["deltae14"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN deltae40
    out["deltae40"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN deltae41
    out["deltae41"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN deltae42
    out["deltae42"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN heating_grainpe
    out["heating_grainpe"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_2level
    out["cooling_2level"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_3level
    out["cooling_3level"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_lya
    out["cooling_lya"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_h2
    out["cooling_h2"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_c0
    out["cooling_c0"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_cplus
    out["cooling_cplus"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_o0
    out["cooling_o0"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_co
    out["cooling_co"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_dust_coll
    out["cooling_dust_coll"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN cooling_dust_rec
    out["cooling_dust_rec"] = $aux_func$
    $$JAFF END
    $$JAFF GET aux_func IN heatingcoolingrate
    out["heatingcoolingrate"] = $aux_func$
    $$JAFF END
    return out
