import sys

# ****************
def parse_prizmo(line):
    """Parse a PRIZMO format reaction line.
    
    Parses a chemical reaction line in PRIZMO format which uses '->' separator
    and optional temperature ranges in square brackets.
    
    Args:
        line (str): A reaction line in PRIZMO format.
                   Example: "H + O2 [10, 1000] -> OH + O 1.2e-10"
    
    Returns:
        tuple: A tuple containing:
            - rr (list): List of reactant species names
            - pp (list): List of product species names  
            - tmin (float): Minimum temperature (K), defaults to 3.0
            - tmax (float): Maximum temperature (K), defaults to 1e6
            - rate (str): Rate coefficient expression as a string
    
    Note:
        - Converts species names: HE->He, E->e-, GRAIN0->GRAIN
        - Replaces variables: user_crflux->crate, user_av->av
        - Adds ',1e10' to PHOTO reactions
    """

    srow = line.strip()

    # temperature ranges
    arow = srow.replace("[", "]").split("]")
    reaction, tlims, rate = [x.strip() for x in arow]

    tmin = tmax = ""
    if "," in tlims:
        tmin, tmax = [x.strip() for x in tlims.split(",")]

    if tmin == "":
        tmin = None
    else:
        tmin = float(tmin.replace("d", "e"))
        if tmin <= 0e0:
            tmin = None

    if tmax == "":
        tmax = None
    else:
        tmax = float(tmax.replace("d", "e"))
        if tmax >= 1e8:
            tmax = None

    reaction = reaction.replace("HE", "He")
    reaction = reaction.replace(" E", " e-")
    reaction = reaction.replace("E ", "e- ")
    reaction = reaction.replace("GRAIN0", "GRAIN")

    rate = rate.replace("user_crflux", "crate")
    rate = rate.replace("user_av", "av")

    rr, pp = reaction.split("->")
    rr = [x.strip() for x in rr.split(" + ")]
    pp = [x.strip() for x in pp.split(" + ")]

    return rr, pp, tmin, tmax, rate

# ****************
def parse_udfa(line):
    """Parse a UDFA format reaction line.
    
    Parses a chemical reaction line in UDFA (UMIST Database for Astrochemistry) format
    which uses ':' as delimiter between fields.
    
    Args:
        line (str): A reaction line in UDFA format with colon-separated fields.
                   Example: "1234:PH:H2:O::H2O::::1.0e-10:0.0:5.0:10.0:1000.0"
    
    Returns:
        tuple: A tuple containing:
            - rr (list): List of reactant species names
            - pp (list): List of product species names
            - tmin (float): Minimum temperature (K)
            - tmax (float): Maximum temperature (K)
            - rate (str): Rate coefficient expression as a string
    
    Note:
        - Reaction types: CR (cosmic ray), PH (photo), standard
        - Skips species: CR, CRP, PHOTON, CRPHOT, empty strings
        - Rate formula depends on reaction type and coefficients ka, kb, kc
    """

    arow = line.split(":")
    rtype = arow[1]
    rr = arow[2:4]
    pp = arow[4:8]
    ka, kb, kc = [float(x) for x in arow[9:12]]
    tmin, tmax = [float(x) for x in arow[12:14]]

    if tmin <= 0e0:
        tmin = None
    if tmax >= 41000.:
        tmax = None

    rate = None
    if rtype == "CR":
        rate = "%.2e * crate" % kc
    elif rtype == "PH":
        rate = "%.2e * exp(-%.2f * av)" % (ka, kc)
    else:
        rate = "%.2e" % ka
        if kb != 0e0:
            rate += " * (tgas / 3e2)**(%.2f)" % kb
        if kc != 0e0:
            rate += " * exp(-%.2f / tgas)" % kc

    skip_species = ["CR", "CRP", "PHOTON", "CRPHOT", ""]

    rr = [x.strip() for x in rr if x.strip() not in skip_species]
    pp = [x.strip() for x in pp if x.strip() not in skip_species]

    return rr, pp, tmin, tmax, rate

# ****************
def parse_kida(line):
    """Parse a KIDA format reaction line.
    
    Parses a chemical reaction line in KIDA (KInetic Database for Astrochemistry) format
    which uses fixed-width columns for different fields.
    
    Args:
        line (str): A reaction line in KIDA format with fixed-width columns.
                   The line should contain reactants, products, and rate coefficients
                   at specific column positions.
    
    Returns:
        tuple: A tuple containing:
            - rr (list): List of reactant species names
            - pp (list): List of product species names
            - tmin (float): Minimum temperature (K)
            - tmax (float): Maximum temperature (K) 
            - rate (str): Rate coefficient expression as a string
    
    Note:
        - Ignores species: CR, CRP, Photon
        - Supports 5 different formula types (1-5) for rate calculations
        - Formula 1: cosmic ray, 2: photo, 3: standard Arrhenius
        - Formula 4-5: special ion-neutral reactions
        - Returns rate="0e0" for unimplemented formulas
    """

    ignore = ["CR", "CRP", "Photon"]

    products_pos = 34
    a_pos = 91

    srow = line

    rr = [x for x in srow[:products_pos].split() if x != '+']
    pp = [x for x in srow[products_pos:a_pos].split() if x != '+']
    arow = srow[a_pos:].split()
    
    # Check if we have enough fields
    if len(arow) < 10:
        raise ValueError(f"Invalid KIDA format: insufficient fields in line: {line}")
    
    ka, kb, kc = [float(x) for x in arow[:3]]
    formula = int(arow[9])
    tmin = float(arow[7])
    tmax = float(arow[8])

    if float(tmin) <= 0e0:
        tmin = None
    if float(tmax) >= 9999.:
        tmax = None


    rate = ""

    if formula == 1:
        rate += "%e * crate" % ka
    elif formula == 2:
        rate += "%.2e * exp(-%e*av)" % (ka, kc)
    elif formula == 3:
        rate += "%.2e" % ka
        if kb != 0e0:
            rate += " * (tgas / 3e2)**(% .2f)" % kb
        if kc != 0e0:
            rate += " * exp(-% .2f / tgas)" % kc
    elif formula == 4:
        rate += "%.2e" % (ka * kb)
        if kc != 0e0:
            rate += " * (0.62 + 0.4767 * %.2e * sqrt(3e2 / tgas))" % kc
    elif formula == 5:
        rate += "%.2e" % (ka * kb)
        if kc != 0e0:
            rate += " * (1e0 + 0.0967 * %.2e * sqrt(3e2 / tgas + %e * 3e2 / 10.526 / tgas))" % (kc, kc**2)
    else:
        print("WARNING: KIDA formula %d not implemented, rate coefficient set to 0e0" % formula)
        rate = "0e0"
        #sys.exit(1)

    rr = [x.strip() for x in rr if x.strip() not in ignore]
    pp = [x.strip() for x in pp if x.strip() not in ignore]

    return rr, pp, tmin, tmax, rate


# ****************
def parse_krome(line, fmt):
    """Parse a KROME format reaction line.
    
    Parses a chemical reaction line in KROME format which uses comma-separated values
    with a format specification header.
    
    Args:
        line (str): A reaction line in KROME format (comma-separated).
        fmt (str): Format specification string starting with '@format:' 
                   that defines the meaning of each comma-separated field.
                   Example: '@format:idx,r,r,p,p,tmin,tmax,rate'
    
    Returns:
        tuple: A tuple containing:
            - rr (list): List of reactant species names
            - pp (list): List of product species names
            - tmin (float): Minimum temperature (K), defaults to 3.0
            - tmax (float): Maximum temperature (K), defaults to 1e6
            - rate (str): Rate coefficient expression as a string
    
    Note:
        - Format fields: r=reactant, p=product, tmin/tmax=temperature limits, rate=rate expr
        - Converts species: E/e->e-, g->"", HE->He
        - Replaces variables: user_crflux/user_crate->crate, user_av->av
        - Converts Fortran syntax to Python using f90_convert()
    """

    line = line.replace(" ", "")

    afmt = [x.strip() for x in fmt.lower().strip().split(":")[1].split(",")]

    arow = [x.strip() for x in line.strip().split(",")]

    assert len(arow) == len(afmt), "ERROR: KROME format does not match line '%s'" % line

    tmin = tmax = None

    tminmax_reps = {"d": "e",
                    ".le.": "",
                    ".ge.": "",
                    ".lt.": "",
                    ".gt.": "",
                    ">": "",
                    "<": ""}

    rr = []
    pp = []
    for i, x in enumerate(afmt):
        if x == "r":
            rr.append(arow[i])
        elif x == "p":
            pp.append(arow[i])
        elif x == "tmin":
            if arow[i].strip().lower() != "none":
                tmin = arow[i].lower()
                for k, v in tminmax_reps.items():
                    tmin = tmin.replace(k, v)
                tmin = float(tmin)
        elif x == "tmax":
            if arow[i].strip().lower() != "none":
                tmax = arow[i].lower()
                for k, v in tminmax_reps.items():
                    tmax = tmax.replace(k, v)
                tmax = float(tmax)
        elif x == "rate":
            rate = arow[i].strip()
        elif x == "idx":
            pass
        else:
            print("ERROR: unknown KROME format %s in line '%s'" % (x, line))
            sys.exit(1)

    rate_reps = {"user_crflux": "crate",
                 "user_crate": "crate",
                 "user_av": "av"}

    for k, v in rate_reps.items():
        rate = rate.replace(k, v)

    rate = f90_convert(rate)

    if "auto" in rate:
        rate = rate.replace("auto", "PHOTO, 1e99")

    sp_reps = {"E": "e-",
                "e": "e-",
                "g": ""}

    rr = [sp_reps[x] if x in sp_reps else x for x in rr]
    pp = [sp_reps[x] if x in sp_reps else x for x in pp]

    sp_sreps = {"HE": "He"}

    for k, v in sp_sreps.items():
        rr = [x.replace(k, v) for x in rr]
        pp = [x.replace(k, v) for x in pp]

    rr = [x.strip() for x in rr if x.strip() != ""]
    pp = [x.strip() for x in pp if x.strip() != ""]

    return rr, pp, tmin, tmax, rate

# ****************
def parse_uclchem(line):
    """Parse a UCLCHEM format reaction line.
    
    Parses a chemical reaction line in UCLCHEM format which uses comma-separated values
    with specific ordering of reactants, products, and rate coefficients.
    
    Args:
        line (str): A reaction line in UCLCHEM format (comma-separated).
                   Format: Reactant1,Reactant2,Reactant3,Product1,Product2,Product3,Product4,
                          Alpha,Beta,Gamma,T_min,T_max,extrapolate
    
    Returns:
        tuple: A tuple containing:
            - rr (list): List of reactant species names
            - pp (list): List of product species names
            - tmin (float): Minimum temperature (K)
            - tmax (float): Maximum temperature (K)
            - rate (str): Rate coefficient expression (currently always "0e0")
    
    Warning:
        This parser is not fully implemented and always returns rate="0e0".
    
    Note:
        - Converts species prefixes: #->_DUST suffix, @->_BULK suffix
        - Handles special reactions: CRP, CRPHOT, PHOTON, FREEZE
        - Ignores many process keywords like ER, DESOH2, THERM, etc.
        - Species conversions: E-->e-, HE->He, SI->Si, CL->Cl, MG->Mg
    """

    srow = line.strip()

    print("WARNING: UCLCHEM reaction format detected, but not implemented yet. Rate set to 0e0.")

    arow = [x.strip() for x in srow.strip().split(",")]
    # !Reactant 1,Reactant 2,Reactant 3,Product 1,Product 2,Product 3,Product 4,Alpha,Beta,Gamma,T_min,T_max,extrapolate
    rr = [x.strip() for x in arow[:3]]
    pp = [x.strip() for x in arow[3:7]]
    ka = arow[7].strip()
    kb = arow[8].strip()
    kc = arow[9].strip()
    extrapolate = arow[12].strip().lower() == "true"

    if extrapolate:
        tmin = 3e0
        tmax = 1e6
    else:
        tmin = float(arow[10])
        tmax = float(arow[11])

    if ",CRP," in srow:
        rate = "%.2e * crate" % float(ka)
    elif ",CRPHOT," in srow:
        rate = "%.2e * (tgas/3e2)**(%.2f) * crate" % (float(ka), float(kb))
    elif ",PHOTON," in srow:
        rate = "%.2e * fuv * exp(-%.2f * av)" % (float(ka), float(kc))
    elif ",FREEZE," in srow:
        rate = "(1e0 + %.2e * 1.671e-3/tgas/asize)*nuth*sigmah*sqrt(tgas/m)" % float(kb)
    else:
        rate = "0e0"

    # FIXME: this is because the parser needs to be finished
    rate = "0e0"


    def convert(sp):

        if sp.startswith("#"):
            sp = sp[1:] + "_DUST"
        if sp.startswith("@"):
            sp = sp[1:] + "_BULK"
        if sp == "E-":
            sp = "e-"

        reps ={"HE": "He",
               "SI": "Si",
               "CL": "Cl",
               "MG": "Mg"}

        for k, v in reps.items():
            sp = sp.replace(k, v)

        return sp

    ignore = ["CR", "CRP", "CRPHOT", "PHOTON", "NAN", ""]
    ignore += ["ER", "ERDES", "FREEZE", "H2FORM", "BULKSWAP", "DESCR",
               "DESOH2", "DEUVCR", "LH", "LHDES", "SURFSWAP", "THERM"]

    rr = [convert(x.strip()) for x in rr if x.strip() not in ignore]
    pp = [convert(x.strip()) for x in pp if x.strip() not in ignore]

    return rr, pp, tmin, tmax, rate

# ****************
def f90_convert(line):
    """Convert Fortran 90 syntax to Python-compatible expressions.
    
    Converts various Fortran 90 syntax elements in rate expressions to
    Python-compatible format for evaluation.
    
    Args:
        line (str): A string containing Fortran 90 syntax elements.
                   Example: "1.0d-10 * dexp(-100d0/T)"
    
    Returns:
        str: The converted string with Python-compatible syntax.
             Example: "1.0e-10 * exp(-100e0/T)"
    
    Note:
        - Converts dexp() to exp()
        - Removes array slice notation (:)
        - Converts Fortran double precision notation (d) to Python scientific notation (e)
          Example: 1.0d-10 becomes 1.0e-10
    """
    import re
    # dexp -> exp
    line = line.replace("dexp(", "exp(")
    line = line.replace("(:)", "")
    # double precision exponential to standard scientific notation
    fa = re.findall(r"[0-9_.]d[0-9_+-]", line)
    #line = re.sub(r"([0-9_.]+)d([0-9_+-]+)", r"\1e\2", line)
    for a in fa:
        line = line.replace(a, a[0]+"e"+a[2])

    return line
