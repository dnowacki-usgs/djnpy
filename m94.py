def m94(ubr, wr, ucr, zr, phiwc, kN, iverbose):
    """M94 - Grant-Madsen model from Madsen(1994)

    function m =  m94( ubr, wr, ucr, zr, phiwc, kN, iverbose )

    Input:
      ubr = rep. wave-orbital velocity amplitude outside wbl [m/s]
      wr = rep. angular wave frequency = 2pi/T [rad/s]
      ucr = current velocity at height zr [m/s]
      zr = reference height for current velocity [m]
      phiwc = angle between currents and waves at zr (radians)
      kN = bottom roughness height (e.q. Nikuradse k) [m]
      iverbose = switch; when 1, extra output
    Returned in structure m:
      m.ustrc  = current friction velocity         u*c [m/s]
      m.ustrr  = w-c combined friction velocity    u*r [m/s]
      m.ustrwm = wave max. friction velocity      u*wm [m/s]
      m.dwc = wave boundary layer thickness [m]
      m.fwc = wave friction factor [ ]
      m.zoa = apparent bottom roughness [m]


    Chris Sherwood, USGS
    November 2005: Removed when waves == 0
    July 2005: Removed bug found by JCW and RPS
    March 2015: Translated to Python by Dan Nowacki, USGS
    """
    import numpy as np
    import math

    MAXIT = 20
    vk = 0.41
    rmu = np.zeros([MAXIT])
    Cmu = np.zeros([MAXIT])
    fwci = np.zeros([MAXIT])
    dwci = np.zeros([MAXIT])
    ustrwm2 = np.zeros([MAXIT])
    ustrr2 = np.zeros([MAXIT])
    ustrci = np.zeros([MAXIT])

    m = {}

    #     ...junk return values
    m["ustrc"] = 99.99
    m["ustrwm"] = 99.99
    m["ustrr"] = 99.99
    fwc = 0.4
    m["fwc"] = fwc
    zoa = kN / 30.0
    m["zoa"] = zoa
    m["dwc"] = kN

    #     ...some data checks
    if wr <= 0.0:
        print("WARNING: Bad value for frequency in M94: wr=", wr)
        return
    if ubr < 0.0:
        print("WARNING: Bad value for orbital vel. in M94: ub=", ubr)
        return
    if kN < 0.0:
        print("WARNING: Wierd value for roughness in M94: kN=", kN)
        return
    if ((zr < zoa) | (zr < 0.05)) & (iverbose == 1):
        print("WARNING: Low value for ref. level in M94: zr=", zr)
        return

    zo = kN / 30.0

    if ubr <= 0.01:
        if ucr <= 0.01:
            #  ...no waves or currents
            ustrc = 0.0
            ustrwm = 0.0
            ustrr = 0.0
            m["ustrc"] = ustrc
            m["ustrr"] = ustrr
            m["ustrwm"] = ustrwm
            m["dwc"] = m["dwc"]
            m["fwc"] = fwc
            m["zoa"] = zoa
            return m

        #  ...no waves
        ustrc = ucr * vk / math.log(zr / zo)
        ustrwm = 0.0
        ustrr = ustrc
        m["ustrc"] = ustrc
        m["ustrr"] = ustrr
        m["ustrwm"] = ustrwm
        m["dwc"] = m["dwc"]
        m["fwc"] = fwc
        m["zoa"] = zoa
        return m

    cosphiwc = abs(math.cos(phiwc))
    rmu[0] = 0.0
    Cmu[0] = 1.0
    fwci[0] = fwc94(Cmu[0], Cmu[0] * ubr / (kN * wr))  # Eqn. 32 or 33
    ustrwm2[0] = 0.5 * fwci[0] * ubr * ubr  # Eqn. 29
    ustrr2[0] = Cmu[0] * ustrwm2[0]  # Eqn. 26
    ustrr = math.sqrt(ustrr2[0])
    dwci[0] = kN
    if (Cmu[0] * ubr / (kN * wr)) >= 8.0:
        dwci[0] = 2 * vk * ustrr / wr
    lnzr = math.log(zr / dwci[0])
    lndw = math.log(dwci[0] / zo)
    lnln = lnzr / lndw
    bigsqr = -1.0 + math.sqrt(1 + ((4.0 * vk * lndw) / (lnzr * lnzr)) * ucr / ustrr)
    ustrci[0] = 0.5 * ustrr * lnln * bigsqr
    nit = 1

    for i in np.arange(1, MAXIT):
        rmu[i] = ustrci[i - 1] * ustrci[i - 1] / ustrwm2[i - 1]
        Cmu[i] = (1 + 2 * rmu[i] * cosphiwc + rmu[i] * rmu[i]) ** 0.5  # Eqn 27
        fwci[i] = fwc94(Cmu[i], (Cmu[i] * ubr / (kN * wr)))  # Eqn. 32 or 33
        ustrwm2[i] = 0.5 * fwci[i] * ubr * ubr  # Eqn. 29
        ustrr2[i] = Cmu[i] * ustrwm2[i]  # Eqn. 26
        ustrr = math.sqrt(ustrr2[i])
        dwci[i] = kN
        if (Cmu[i] * ubr / (kN * wr)) >= 8.0:
            dwci[i] = 2 * vk * ustrr / wr  # Eqn.36
        lnzr = math.log(zr / dwci[i])
        lndw = math.log(dwci[i] / zo)
        lnln = lnzr / lndw
        bigsqr = -1 + math.sqrt(1 + ((4.0 * vk * lndw) / (lnzr * lnzr)) * ucr / ustrr)
        ustrci[i] = 0.5 * ustrr * lnln * bigsqr
        # Eqn. 38
        diffw = abs((fwci[i] - fwci[i - 1]) / fwci[i])
        if diffw < 0.0005:
            break
        nit = nit + 1

    ustrwm = math.sqrt(ustrwm2[nit])
    ustrc = ustrci[nit]
    ustrr = math.sqrt(ustrr2[nit])

    zoa = math.exp(
        math.log(dwci[nit]) - (ustrc / ustrr) * math.log(dwci[nit] / zo)
    )  # Eqn. 11
    fwc = fwci[nit]
    dwc = dwci[nit]

    if iverbose == 1:
        for i in range(nit):
            print(
                "i=",
                i,
                "fwc=",
                fwci[i],
                "dwc=",
                dwci[i],
                "u*c=",
                ustrci[i],
                "u*r=",
                math.sqrt(ustrwm2[i]),
                "u*r=",
                math.sqrt(ustrr2[i]),
            )

    m["ustrc"] = ustrc
    m["ustrr"] = ustrr
    m["ustrwm"] = ustrwm
    m["dwc"] = dwc
    m["fwc"] = fwc
    m["zoa"] = zoa

    return m


def fwc94(cmu, cukw):
    # FWC94 - Wave-current friction factor
    # Equations 32 and 33 in Madsen, 1994

    import math

    fwc = 0.00999  # meaningless (small) return value

    if cukw <= 0.0:
        print("ERROR: cukw too small in fwc94:", cukw)
        return

    if cukw < 0.2:
        fwc = math.exp(7.02 * 0.2 ** (-0.078) - 8.82)
        print("WARNING: cukw very small in fwc94:", cukw)

    if (cukw >= 0.2) & (cukw <= 100.0):
        fwc = cmu * math.exp(7.02 * cukw ** (-0.078) - 8.82)
    elif (cukw > 100.0) & (cukw <= 10000.0):
        fwc = cmu * math.exp(5.61 * cukw ** (-0.109) - 7.30)
    elif cukw > 10000:
        fwc = cmu * math.exp(5.61 * 10000 ** (-0.109) - 7.30)
    else:
        print("WARNING: cukw very large in fwc94:", cukw)

    return fwc


# pi = 3.14159
# m = m94(0.1, 2*pi/5, .3, 1, 0, 0.05, 1)
# print m
