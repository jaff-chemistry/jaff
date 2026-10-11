# ABOUTME: jaffgen test template: every auxiliary function of the GOW network
# ABOUTME: Rendered via GET aux_func; compared numerically against tests/golden
import math


def aux_functions():
    out = {}
    # $JAFF GET aux_func FOR kcr_h_fac
    out["kcr_h_fac"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR qcr_h
    out["qcr_h"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR qcr_h2
    out["qcr_h2"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR qcr
    out["qcr"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR ncrh2
    out["ncrh2"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR kgr_gong
    out["kgr_gong"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR n2ncr
    out["n2ncr"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR fshield_h2
    out["fshield_h2"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR fshield_c
    out["fshield_c"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate0
    out["chemrate0"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate1
    out["chemrate1"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate4
    out["chemrate4"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate8
    out["chemrate8"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate10
    out["chemrate10"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate13
    out["chemrate13"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate14
    out["chemrate14"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate15
    out["chemrate15"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate16
    out["chemrate16"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate17
    out["chemrate17"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate18
    out["chemrate18"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate19
    out["chemrate19"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate20
    out["chemrate20"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate21
    out["chemrate21"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate22
    out["chemrate22"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate23
    out["chemrate23"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate27
    out["chemrate27"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate28
    out["chemrate28"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate32
    out["chemrate32"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate35
    out["chemrate35"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate37
    out["chemrate37"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate39
    out["chemrate39"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate40
    out["chemrate40"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate41
    out["chemrate41"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate42
    out["chemrate42"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate48
    out["chemrate48"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR chemrate49
    out["chemrate49"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR deltae0
    out["deltae0"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR deltae1
    out["deltae1"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR deltae2
    out["deltae2"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR deltae13
    out["deltae13"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR deltae14
    out["deltae14"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR deltae40
    out["deltae40"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR deltae41
    out["deltae41"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR deltae42
    out["deltae42"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR heating_grainpe
    out["heating_grainpe"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_2level
    out["cooling_2level"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_3level
    out["cooling_3level"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_lya
    out["cooling_lya"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_h2
    out["cooling_h2"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_c0
    out["cooling_c0"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_cplus
    out["cooling_cplus"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_o0
    out["cooling_o0"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_co
    out["cooling_co"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_dust_coll
    out["cooling_dust_coll"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR cooling_dust_rec
    out["cooling_dust_rec"] = $aux_func$
    # $JAFF END
    # $JAFF GET aux_func FOR heatingcoolingrate
    out["heatingcoolingrate"] = $aux_func$
    # $JAFF END
    return out
