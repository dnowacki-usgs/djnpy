from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity


@ureg.wraps(
    ureg.gram / ureg.cm**3,
    (ureg.gram, ureg.gram, ureg.gram / ureg.cm**3, ureg.gram / ureg.cm**3, None),
)
def drybulkdens(
    water_mass,
    unadjusted_dry_sediment_mass,
    water_density=Q_(1, "g/cm**3"),
    sediment_density=Q_(2.65, "g/cm**3"),
    salinity=35,
):
    """Compute dry bulk density

    Parameters
    ----------
    water_mass :
        Water mass, with units as specified. If no units specified, defaults to gram / centimeter ** 3
    unadjusted_dry_sediment_mass
        Dry sediment mass (unadjusted for dried salts), with units as specified. If no units specified, defaults to gram / centimeter ** 3
    water_density : pint.util.Quantity, default 1 gram / centimeter ** 3
        Water density
    sediment_density : pint.util.Quantity, default 2.65 gram / centimeter ** 3
        Sediment density
    salinity : float, default 35
        Salinity

    Returns
    -------
    pint.util.Quantity
        Dry bulk density in gram / centimeter ** 3
    """

    dry_sediment_mass = adjust_mass_for_salt(
        water_mass, salinity, unadjusted_dry_sediment_mass
    )

    n = porosity(water_mass, water_density, dry_sediment_mass, sediment_density)

    return (1 - n) * sediment_density


@ureg.wraps(
    ureg.gram / ureg.cm**3,
    (ureg.gram, ureg.gram, ureg.gram / ureg.cm**3, ureg.gram / ureg.cm**3, None),
)
def wetbulkdens(
    water_mass,
    unadjusted_dry_sediment_mass,
    water_density=Q_(1, "g/cm**3"),
    sediment_density=Q_(2.65, "g/cm**3"),
    salinity=35,
):
    """Compute wet bulk density

    Parameters
    ----------
    water_mass :
        Water mass, with units as specified. If no units specified, defaults to gram / centimeter ** 3
    unadjusted_dry_sediment_mass
        Dry sediment mass (unadjusted for dried salts), with units as specified. If no units specified, defaults to gram / centimeter ** 3
    water_density : pint.util.Quantity, default 1 gram / centimeter ** 3
        Water density
    sediment_density : pint.util.Quantity, default 2.65 gram / centimeter ** 3
        Sediment density
    salinity : float, default 35
        Salinity

    Returns
    -------
    pint.util.Quantity
        Wet bulk density in gram / centimeter ** 3
    """

    dry_sediment_mass = adjust_mass_for_salt(
        water_mass, salinity, unadjusted_dry_sediment_mass
    )

    return (dry_sediment_mass + water_mass) / (
        dry_sediment_mass / sediment_density + water_mass / water_density
    )


def adjust_mass_for_salt(water_mass, salinity, unadjusted_dry_sediment_mass):

    salt_mass = water_mass * salinity / 1000
    return unadjusted_dry_sediment_mass - salt_mass


def porosity(water_mass, water_density, dry_sediment_mass, sediment_density):

    return (water_mass / water_density) / (
        dry_sediment_mass / sediment_density + water_mass / water_density
    )


# print(Q_(1, "g/cm**3"))
# print(wetbulkdens(1 * ureg.gram, 2 * ureg.gram, salinity=20))
# print(drybulkdens(1, 2 * ureg.gram, salinity=20))
