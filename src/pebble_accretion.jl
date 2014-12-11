module PA

using MKS
# using Roots
using Gadfly

# Set global parameters
ALPHA_VISCOSITY = 1e-4
PEBBLE_SIZE = 1e0
MIN_PROTOPLANET_SIZE = 1e3
MAX_PROTOPLANET_SIZE = 2e3
PROTOPLANET_DENSITY = 3e0
PEBBLE_DENSITY = 2e0
β0 = 5e2 #Initial gas surface density at 1AU [g/cm^2]
τDISSIPATION = 3e6  #The disk gas dissipation timescale
DUSTTOGAS = 1e-2  #Initial dust to gas ratio
INITIAL_SURFACE_DENSITY = 1  # Initial surface density of planetesimals+protoplanets at 1AU
ALPHA_SLOPE = -1.5  # Initial slope of the surface density
PLANETESIMAL_SIZE = 1e2



############################################################################
# General Functions
############################################################################

function set_globals(;pebble_density = 1e0,
                     α_t = 1e-4,
                     pebble_size = 1e0,
                     min_protoplanet_size = 1e3,
                     max_protoplanet_size = 2e3,
                     protoplanet_density = 3e0,
                     β_0 = 5e2, #Initial gas surface density at 1AU [g/cm^2]
                     τ_dissipation = 3e6,  #The disk gas dissipation timescale
                     dust_to_gas = 1e-2,  #Initial dust to gas ratio
                     initial_surface_density = 1e0,  # Initial surface density of planetesimals+protoplanets at 1AU
                     α_slope = -1.5, # Initial slope of the surface density)
                     planetesimal_size =1e2)
    global PEBBLE_DENSITY = pebble_density
    global ALPHA_VISCOSITY = α_t
    global PEBBLE_SIZE = pebble_size
    global MIN_PROTOPLANET_SIZE = min_protoplanet_size
    global MAX_PROTOPLANET_SIZE = max_protoplanet_size
    global PROTOPLANET_DENSITY = protoplanet_density
    global β0 = β_0
    global τDISSIPATION = τ_dissipation
    global DUSTTOGAS = dust_to_gas
    global INITIAL_SURFACE_DENSITY = initial_surface_density
    global ALPHA_SLOPE = α_slope
    global PLANETESIMAL_SIZE = planetesimal_size
end

function stokes(particle_size, a, t; ρ=2e0)
    """Calculate the stokes number.

        This calculation assumes models for the gas density, sound speed and
        mean free path:
            ρ_g = Σ / H   where H/r = 0.033 * a^1/4
            Σ=β_0 * exp(t/3e6years) / a
            c_s = 0.033 * a^1/4 * v_k
            λ = 2e-9 / ρ_g (cgs units)

        Args:
            particle_size: Size of accreting particles (cm).
            a: Semi-major axis (AU).
            t: Time in years.
            ρ: Material density of particle.
    """

    #First calculate the stopping time
    ρ_g = gasdensity(a, t)
    v_k = (GRAV*MSUN/a/AU)^0.5
    c = 0.033 * a^0.25 * v_k * 100.  # sound speed converted to cgs
    λ = 2e-9 / ρ_g  # cgs units only (from OK10)
    s_max = 9000. * c / 1e5 * 1e-10 / ρ_g * 30 / headwindvelocity(a) # cm
    s_transition = 9. * λ / 4.
    if particle_size < s_transition  # Epstein drag
        # println("Epstein")
        t_s = ρ * particle_size / ρ_g / c
    elseif particle_size < s_max   # Stokes drag
        # println("Stokes")
        t_s = (4ρ * particle_size * particle_size /
               9ρ_g / c / λ)
    else  # Quadratic drag
        # println("Quadratic")
        t_s = 6ρ*particle_size/ρ_g/headwindvelocity(a)/100
        # error("Particle too large: Calculations not valid")
    end

    #Now calculate stokes number
    omega_0 = v_k / a / AU
    t_s * omega_0
end


function hill_radius(a, r, ρ)
    """ Calculate the hill radius.

        Args:
            a: Semi-major axis in AU.
            r: Planet radius in km.
            ρ: Planet density in g/cm^3
        Return:
            Hill radius in meters.

    """

    mass = 4./3. * π * (r*1000.*100)^3 * ρ / 1000 #kg
    a*(mass/3MSUN)^(1./3.)*AU
end


function sphere_mass(r, ρ)
    """ Calculate the mass.

        Args:
            r: Planet radius in km.
            ρ: Planet density in g/cm^3
        Return:
            mass in kg.

    """

    4./3. * π * (r*1000.*100)^3 * ρ / 1000 #kg
end


function sphere_radius(mass, ρ)
    """ Calculate the radius of corresponding sphere.

        Args:
            mass: Mass of object in Earth masses.
            ρ: density of object in g/cm^3.
        Return:
            Object radius in km.

    """
    (mass * MEARTH / GRAM / ρ / 4 * 3 / π)^(1./3) / 100 / 1000  #km
end


function gasdensity(a, t)
    """Calculate gas density as a function of semi-major axis.

        Args:
            a: Semi-major axis (AU).
            t: Time in years.
        Return:
            Gas density in g/cm^3.

    """

    scale_height = 0.033(a^1.25) * AU * a * 100  # cm
    sigma = β0 * exp(-t/τDISSIPATION) / a  # g/cm^2
    sigma / scale_height / 2
end


function gassurfacedensity(a, t)
    """Calculate gas surface density as a function of semi-major axis.

        Args:
            a: Semi-major axis (AU).
            t: Time in years.
        Return:
            Gas density in g/cm^3.

    """

    β0 * exp(-t/τDISSIPATION) / a  # g/cm^2
end


function headwindvelocity(a)
    """ Calculate the headwind velocity of a protplanet.

        Args:
            a: Semi-major axis (AU).
        Return:
            Velocity of headwind in m/s.
    """

    v_k = (GRAV*MSUN/a/AU)^ 0.5
    # c = 0.033 * a^0.25 * v_k  # From LJ12 Eq. 3
    # c * c / v_k
    0.0015 * sqrt(a) * v_k # LJ14 eq 18
end


function eta(a)
    """ Measure of radial gas support. From eqn 18 of LJ14.
        Args:
            a: Semi-major axis in AU.
        Return:
            eta: dimensionless
    """

    0.0015*sqrt(a)
end


function mass2surfacedensity(mass, inner, outer)
    """ Calculate the surface density in an annulus.

        Args:
            mass: Mass in annulus [M_Earth].
            inner: Inner edge of annulus in AU.
            outer: Outer edge of annulus in AU.

        Return:
            Surface desnity in annulus in g/cm^2.

    """


    area = π*((outer*100AU)^2 - (inner*100AU)^2)
    mass * 1000MEARTH / area
end


function surfacedensity2mass(Σ, inner, outer; α=0)
    """ Calculate the mass in an annulus.

        Args:
            Σ: Surface density in annulus [g/cm^2].
            inner: Inner edge of annulus in AU.
            outer: Outer edge of annulus in AU.
            α: Slope of surface density. (Σ = Σ_0 * (a/AU)^α)

        Return:
            Mass in annulus in M_Earth.

    """
    Σ_mearth_per_au2 = Σ/1000MEARTH*100AU*100AU
    2π*Σ_mearth_per_au2/(α+2)*(outer^(α+2)-inner^(α+2))
end


############################################################################
# Pebble Production and Surface Density
############################################################################

function pebbleproduction(t; β=β0, m_star=1., Z=DUSTTOGAS)
    """ Determine the rate at which migrating pebbles are produced in the outer disk.
       (Eqn. 14 from Lambrechts and Johansen 2014)

        Args:
            t: Time in years.
            β: Gas surface density at 1 AU (g/cm^2).
            m_star: Mass of central star (M_sun).
            Z: Initial dust to gas ratio.
        Return:
            Pebble productino rate. (M_earth/year)
    """
    # println("here ", t)


    if t < 1e4
        println("Timescales ≲ 1e4 years are not relevant. Setting pebble production disk age to 1e4 years")
        println(t)
        t = 1e4
    end
    9.5e-5(β*exp(-t/τDISSIPATION)/500.) * m_star^(1./3.) * (Z/0.01)^(5./3.) * (t/1e6)^(-1./3.)
end


function pebbledensity(a, t; r0=PEBBLE_SIZE, ρ=PEBBLE_DENSITY, pebble_flux=nothing)
    """ Determine the surface density of pebbles throughout disk.

        Args:
            t: time in years.
            a: semi-major axis in AU.
            r0: dominant pebble size in cm.
            ρ: material density of pebbles
            pebble_flux: The rate at which pebbles enter the annulus
                at a. This allows for a radial change in the pebble
                flux in the disk. [M_earth/yr]

        Return:
            Surface density of pebbles in g/cm^2

    """
    ρ = PEBBLE_DENSITY
    if pebble_flux == nothing
        pebble_flux = pebbleproduction(t)
    end
    r0 = min(dominantpebblesize(a, t), r0)

    pebble_flux * 1000MEARTH/YEAR / (4π*a*100AU*stokes(r0, a, t, ρ=ρ)*100headwindvelocity(a))
    # pebble_flux * 1000MEARTH/YEAR / (4π*a*100AU*st*100headwindvelocity(a))
end


function pebbledensity_LJ(a, t)
    """ Pebble surface density according to eqn 25 in LJ2014

        I don't know what is going on in eqn 25. The two expressions don't seem to be
        that close. (factor of a few)

        Args:
            t: time in years.
            a: semi-major axis in AU.
        Return:
            Pebble surface density in g/cm^2

    """

    m_star = 1
    β = β0 * exp(-t/τDISSIPATION)
    z_0 = 0.01
    # println(0.069 * (β/500) * (z_0/0.01)^(5./6.) * m_star^(-1./12.) * (t/1e6)^(-1./6.) * (a/10.)^(-3./4))

    v_k = (GRAV*MSUN/a/AU)^0.5
    Ω = v_k / a / AU
    Σ = gassurfacedensity(a, t)
    2.^(5./6.) * 3.^(-7./12.) * 0.5^(1./3.) / 0.5^.5 * 0.01^(5./6.) * Σ * (Ω*t*YEAR)^(-1./6.)
    0.069 * (β/500) * (z_0/0.01)^(5./6.) * m_star^(-1./12.) * (t/1e6)^(-1./6.) * (a/10.)^(-3./4)

    sqrt(2*pebbleproduction(t)*1000MEARTH/YEAR * Σ/(sqrt(3.)*π*0.5*100AU*a*100v_k))
end


function dominantpebblesize(a, t)
    """ Dominant pebble size according to eqn 20 of LJ14.

        Args:
            a: semi-major axis in AU.
            t: time in years.
        Return:
            pebble size in cm.

    """

    v_k = (GRAV*MSUN/a/AU)^0.5
    ρ = PEBBLE_DENSITY
    ρ_g = gasdensity(a, t)
    c = 0.033 * a^0.25 * v_k * 100.  # sound speed converted to cgs
    Ω = v_k / a / AU
    η = 0.0015*sqrt(a)
    st = sqrt(3)/16./η*pebbledensity_LJ(a,t)/gassurfacedensity(a,t)
    size = st/Ω/ρ*ρ_g*c
    # st
end


######################################################################################
# Collisional Cross Sections (OK10)
######################################################################################

function alpha_p(a, ρ_s=PROTOPLANET_DENSITY):
    """Calcualate the dimensionless planet size (planet radius over hill radius).

        Args:
            a: Semi-major axis (AU).
            rho_s: Planet density (g/cm^3)

    """
    5.7e-3*(ρ_s/3)^(-1./3.)/a
end


function zeta_w(rho_s, R_p, a):
    """Calculate the dimensionless parameter for headwind velocity.

        Args:
            rho_s: Density of protoplanet (gm/cm^3)
            R_p: Radius of protoplanet (km)
            a: Semi-major axis (AU)
        Return:
            Dimensionless parameter zeta.

    """

    v_hw = headwindvelocity(a)
    12.5 * rho_s^(-1./3.) * v_hw / 30. * 100. / R_p * sqrt(a)
end


function positive_zero(f)
    tol = 1e-10
    error = 1e100
    best_min = 0
    best_max = 1e5
    x=0
    while abs(error) > tol
        x = (best_max + best_min)/2
        error = f(x)
        if error > 0
            best_max = x
        else
            best_min = x
        end
    end
    x
end


function b_set(ζ, st):
    """Calcualate impact parameter in settling regime.
        args:
            ζ: Dimensionless headwind velocity.
            st: Stokes number.

    """

    f(x) = x^3 + x^2*(2.*ζ/3.) - (8*st)
    # roots = fzeros(f)
    # b = maximum(roots)
    b = positive_zero(f)
    # println(b)
    # b=.1
    st_crit = 12. / ζ^3
    b*exp(-((st/st_crit)^0.65))
    # f
end


function b_hyp(ζ, st, α):
    """Calcualate impact parameter in hyperbolic regime.
        args:
            ζ: Dimensionless headwind velocity.
            st: Stokes number.
            α: Dimensionless planet size.

    """

    v_a = ζ * (1.+4.*st*st)^0.5 / (1.+st*st)
    α * (1.+6./α/v_a/v_a)^0.5
end


function b_3b(α, st):
    """Calcualate impact parameter in 3-body regime.
        args:
            st: Stokes number.
            α: Dimensionless planet size.

    """

    1.7 * sqrt(α) + 1.0 / st
end


function impact_parameter_nondimensional(α, ζ, st):
    """Determine regime and calculate appropriate imact parameter.

    Args:
        ζ: Dimensionless headwind velocity.
        st: Stokes number.
        α: Dimensionless planet size.

    """
    if st < min(1., 12./ζ/ζ/ζ)
        # println("settling regime")
        b = maximum([b_set(ζ, st), α])
    elseif st > maximum([ζ, 1.])
        # println("Three body regime")
        b = maximum([b_3b(α, st), α])
    else
        # println("Hyperbolic regime")
        b = maximum([b_set(ζ, st), b_hyp(ζ, st, α)])
    end
    b
end


function accretion_regime_nondimensional(α, ζ, st)
    """ Determine the acccretion regime for the given parameters.

        Args:
        ζ: Dimensionless headwind velocity.
        st: Stokes number.
        α: Dimensionless planet size.

    """
    if st < min(1., 12./ζ/ζ/ζ)
        b = "settling"
    elseif st > maximum([ζ, 1.])
        b = "three body"
    else
        b = "hyperbolic"
    end
    b
end


function accretion_regime(a, t, r_p, r_0; ρ_pl=PROTOPLANET_DENSITY, ρ_peb=PEBBLE_DENSITY)
    """ determine accretion regime from physical units
        Args:
            a: Semi-major axis in AU.
            t: Time in years.
            r_p: Protoplanet radius in km.
            r_0: Pebble size in cm.
            ρ_pl: Planet density in g/cm^3.
            ρ_peb: Pebble density in g/cm^3.
        Return:
            Accretion regime (string).

    """
    α = alpha_p(a, ρ_pl)
    ζ = zeta_w(ρ_pl, r_p, a)
    st = stokes(r_0, a, t, ρ=ρ_peb)
    accretion_regime_nondimensional(α, ζ, st)
end


function impact_parameter(a::FloatingPoint, t::FloatingPoint, r_p::FloatingPoint, r_0::FloatingPoint; ρ_pl=PROTOPLANET_DENSITY, ρ_peb=PEBBLE_DENSITY)
    """ Calculate the impact paramter from and in physical units.

        Args:
            a: Semi-major axis in AU.
            t: Time in years.
            r_p: Protoplanet radius in km.
            r_0: Pebble size in cm.
            ρ_pl: Planet density in g/cm^3.
            ρ_peb: Pebble density in g/cm^3.
        Return:
            Impact parameter in meters.
    """

    α = alpha_p(a, ρ_pl)
    ζ = zeta_w(ρ_pl, r_p, a)
    st = stokes(r_0, a, t, ρ=ρ_peb)
    b = impact_parameter_nondimensional(α, ζ, st)
    r_hill = hill_radius(a, r_p, ρ_pl)
    ans = b*r_hill
    ans
end


########################################################################
# Accretion rate
########################################################################

function relative_velocity(a, b, r0, t; α_t=ALPHA_VISCOSITY, e=0)
    """Calculate the relative velocity between pebbles and protoplanet.

        Args:
            a: Semi-major axis in AU.
            b: Impact parameter in meters.
            r0: Pebble size in cm.
            t: Time in years.
            α_t: Alpha visocity parameter.
            e: Eccentricity.
        Return:
            Relative velocity in m/s.
    """

    v_k = (GRAV*MSUN/a/AU)^0.5
    st = stokes(r0, a, t)
    c = 0.033 * a^0.25 * v_k
    v_turb = sqrt(α_t * c *c * minimum([2st, 1./(1.+st)]))
    v_drift = eta(a)*v_k
    v_ecc = e*v_k
    v_shear = b*v_k/a/AU
    maximum([v_turb, v_drift, v_ecc, v_shear])
end


function growth_rate(a, t, r_p, r0; e=0., α_t=ALPHA_VISCOSITY, ρ_pl=PROTOPLANET_DENSITY, ρ_peb=PEBBLE_DENSITY, pebble_flux=nothing)
    """ Calculate growth rate of a protoplanet.

        Args:
            a: Semi-major axis in AU.
            t: Time in years.
            r_p: Radius of protoplanet in km.
            r0: Size of pebbles in cm.
            e: Eccentricity.
            α_t: Viscocity parameter.
            ρ_pl: Material density of protoplanet in g/cm^3.
            ρ_peb: Material density of pebbles in g/cm^3.
            pebble_flux: Flux of pebbles through disk at a. (Defaults
                to the production rate.)
        Return:
            Growth rate of protoplanet in Earth masses per year.

    """
    b = impact_parameter(a, t, r_p, r0, ρ_pl=ρ_pl, ρ_peb=ρ_peb)

    v_rel = relative_velocity(a, b, r0, t, α_t=α_t, e=e)
    h_gas = 0.033(a^1.25) * AU * a  # meters
    h_solid = h_gas*sqrt(α_t/(α_t+stokes(r0, a, t)))
    # println(h_solid/1000.)
    Σ = pebbledensity(a, t, r0=r0, ρ=ρ_peb, pebble_flux=pebble_flux) * GRAM / CM^2 # MKS
    # Σ = pebbledensity_LJ(a, t)
    if h_solid > b
        # println("3-D")
        ρ = Σ / 2h_solid
        Ṁ = π*b*b*ρ*v_rel * YEAR / MEARTH
    else
        # println("2-D")
        Ṁ = 2*b*Σ*v_rel * YEAR / MEARTH
    end
    # println("Accretion efficiency: ", Ṁ/pebbleproduction(t))
    Ṁ
end


###################################################################
# Planetesimal mass function
###################################################################

function total_accretion_rate(a, t, total_mass; num_bodies=0, body_size=0, pebble_flux=nothing)
    """ Calculate the total accretion rate of pebbles onto all num_bodies bodies.

        This function assumes a few paraeters but is just for getting a sense of how
            the mass distribution of planetesimals affects their total accretion rate.
            A few notes: the total accretion rate should be much less than the pebble
            production rate otherwise this is not valid. The total mass should also be
            low enough to correspond to a relatively small annulus.

        Args:
            a: Semi-major axis in AU.
            t: Time in years.
            total_mass: Total mass of planetesimals in Earth masses.
            num_bodies: Number of planetesimals.
            body_size: Size of planetesimals in km.

        Return:
            The total accretion rate onto all planetesimals in MEARTH/YEAR

    """
    if num_bodies==0 && body_size==0
        error("Must specify number of bodies or body size")
    end
    if num_bodies==0
        individual_mass = sphere_mass(body_size, PROTOPLANET_DENSITY) / MEARTH
        num_bodies = total_mass/individual_mass
        # println(num_bodies)
    else
        individual_mass = total_mass/num_bodies
        body_size = sphere_radius(individual_mass, PROTOPLANET_DENSITY)
    end
    Ṁ_individual = growth_rate(a, t, body_size, PEBBLE_SIZE, pebble_flux=pebble_flux)
    # println(num_bodies, " ", individual_mass/Ṁ_individual, " ", pebble_flux)
    Ṁ_individual * num_bodies
end


#######################################################################
# Full Simulation
#######################################################################

abstract Accreter


type PlanetesimalFlock <: Accreter
    """ Accreter that represents an annulus of planetesimals.

        Feilds:
            sma: Inner extent of annulus [AU].
            outer_extent: Outer extent of annulus [AU].
            mass: Total mass of planetesimals in annulus [M_Earth].
            body_size: Size of planetesimals [km].

    """
    sma::FloatingPoint
    outer_extent::FloatingPoint
    mass::FloatingPoint
    body_size::FloatingPoint
end


type Protoplanet <: Accreter
    """ Accreter that represents a single protoplanet.

        Feilds:
            sma: Semi-major axis of protoplanet [AU].
            mass: Mass of protoplanet [M_Earth].
            eccentricity: Eccentricity of protoplanet.

    """
    sma::FloatingPoint
    mass::FloatingPoint
    eccentricity::FloatingPoint
end


type Disk
    """ Container for all the accreting bodies in the disk.

        Feilds:
            max_surface_density: The maximum surface density of planetesimals before
                a new protoplanet is formed from them [g/cm^2].
            bodies: All the accreting bodies.
            planetesimals: Only the PlanetesimalFlock bodies.
            protoplanets: Only the Protoplanet bodies.

    """
    max_surface_density::FloatingPoint
    bodies::Array{Accreter}
    planetesimals::Array{PlanetesimalFlock}
    protoplanets::Array{Protoplanet}
end




function Disk(inner_edge::FloatingPoint, outer_edge::FloatingPoint, num_annuli::Integer, max_surface_density)
    """ Constructor for type PlanetesimalDisk.

        Args:
            inner_edge: The inner edge of the planetesimals disk in AU.
            outer_edge: The outer edge of the planetesimals disk in AU.
            num_annuli: Defines resolution of simulation.
            max_surface_density: The threshold at which additional mass
                is converted to a new planetary embryo.

        Returns:
            A new PlanetesimalDisk object.

    """

    #Starting with a surface density of Σ = .1*(a/AU)^-1.5
    α = ALPHA_SLOPE
    planetesimals = PlanetesimalFlock[]
    protoplanets = Protoplanet[]

    annuli = [10^loga for loga in linspace(log10(inner_edge), log10(outer_edge), num_annuli+1)]

    # Add the planetesimals
    for i = 1:num_annuli
        mass_in_annulus = surfacedensity2mass(INITIAL_SURFACE_DENSITY/2., annuli[i], annuli[i+1], α=α)
        planetesimal_ring = PlanetesimalFlock(annuli[i], annuli[i+1], mass_in_annulus, PLANETESIMAL_SIZE)
        append!(planetesimals, [planetesimal_ring])
    end

    # Add the protoplanets
    sma = inner_edge
    i=0
    while sma < outer_edge
        radius = rand()*(MAX_PROTOPLANET_SIZE - MIN_PROTOPLANET_SIZE) + MIN_PROTOPLANET_SIZE
        mass = sphere_mass(radius, PROTOPLANET_DENSITY)/MEARTH
        Σ_mearth_per_au2 = INITIAL_SURFACE_DENSITY/2./1000MEARTH*100AU*100AU
        next_sma = (mass*(α+2)/2π/Σ_mearth_per_au2 + sma^(α+2))^(1./(α+2))
        append!(protoplanets, [Protoplanet(rand()*(next_sma - sma) + sma, mass, 0)])
        sma = next_sma
        i+=1
    end

    Disk(max_surface_density, [planetesimals, protoplanets], planetesimals, protoplanets) #random embryos for now
end

import Base.sort!
function sort!(disk::Disk; by=nothing)
    sortkey(accreter) = accreter.sma
    sort!(disk.bodies, by=sortkey)
    sort!(disk.planetesimals, by=sortkey)
    sort!(disk.protoplanets, by=sortkey)
end



# function Disk(inner_edge::FloatingPoint, outer_edge::FloatingPoint,
#               total_mass::FloatingPoint, min_sep::FloatingPoint,
#               max_sep::FloatingPoint)
#     """ Constructs an EmbryoDisk with the provided parameters.

#         Args:
#             inner_edge: Inner extent of disk in AU.
#             outer_edge: Outer extent of disk in AU.
#             total_mass: Total mass of embryos in Earth masses.
#             min_sep: Minimum separation between embryos in mutual hill radii.
#             max_sep: Maximum separation between embryos in mutual hill radii.

#     """
#     x = "I'll leave this one for later"
# end


function accrete_pebbles!(disk::Disk, time, dt)
    """ Calculate the pebble accretion rate at each annulus and apply it for the provided time interval.

        Args:
            disk: <Disk> The protoplanetary disk for which to calculate the pebble accretion rates.
            time: The time in years.
            dt: The time step in years.

    """
    sort!(disk)
    pebble_flux = pebbleproduction(time)
    printed=0
    for i=length(disk.bodies):-1:1
        accretion_rate = accrete_pebbles!(disk.bodies[i], time, dt, pebble_flux)
        pebble_flux -= accretion_rate
    end
    convert_planetesimals_to_embryos!(disk)
end

function accrete_pebbles!(protoplanet::Protoplanet, time, dt, pebble_flux)
    radius = sphere_radius(protoplanet.mass, PROTOPLANET_DENSITY)
    accretion_rate = growth_rate(protoplanet.sma, time, radius, PEBBLE_SIZE, pebble_flux=pebble_flux)
    accretion_rate = min(accretion_rate, pebble_flux)
    protoplanet.mass += accretion_rate * dt
    accretion_rate
end


function accrete_pebbles!(planetesimals::PlanetesimalFlock, time, dt, pebble_flux)
    accretion_rate = total_accretion_rate(planetesimals.sma, time, planetesimals.mass,
                                          body_size=planetesimals.body_size, pebble_flux=pebble_flux)
    accretion_rate = min(accretion_rate, pebble_flux)
    planetesimals.mass += accretion_rate * dt
    # println("planetesimal ", accretion_rate/pebble_flux)
    accretion_rate
end


function convert_planetesimals_to_embryos!(disk::Disk)
    # If there is too much mass in planetesimals create a new embryo
    sort!(disk)
    accumulated_mass = 0
    inner_bound = disk.planetesimals[1].sma
    included_annuli = Integer[]
    excess_Σ = FloatingPoint[]
    for i in 1:length(disk.planetesimals)
        p = disk.planetesimals[i]
        Σ = mass2surfacedensity(p.mass, p.sma, p.outer_extent)
        if Σ > disk.max_surface_density
            append!(included_annuli, [i])
            append!(excess_Σ, [Σ-disk.max_surface_density])
            accumulated_mass += surfacedensity2mass(Σ-disk.max_surface_density, p.sma, p.outer_extent)
            if accumulated_mass > sphere_mass(MIN_PROTOPLANET_SIZE, PROTOPLANET_DENSITY)/MEARTH
                num_bodies = iceil(accumulated_mass/sphere_mass(MAX_PROTOPLANET_SIZE, PROTOPLANET_DENSITY)/MEARTH)
                for j in 1:num_bodies
                    new_embryo = Protoplanet(rand()*(p.outer_extent - inner_bound)+inner_bound, accumulated_mass/num_bodies, 0)
                    append!(disk.protoplanets, [new_embryo])
                    append!(disk.bodies, [new_embryo])
                end
                for j in 1:length(included_annuli)
                    index = included_annuli[j]
                    disk.planetesimals[index].mass -= surfacedensity2mass(excess_Σ[j], disk.planetesimals[index].sma, disk.planetesimals[index].outer_extent)
                end
                included_annuli = Integer[]
                excess_Σ = FloatingPoint[]
                inner_bound = p.outer_extent
            end
        end
    end
end

import Gadfly.plot
function plot(disk::Disk; binned=false, title="")
    sort!(disk)
    planetesimal_sma = Number[]
    planetesimal_Σ = Number[]
    for p in disk.planetesimals
        Σ = mass2surfacedensity(p.mass, p.sma, p.outer_extent)
        append!(planetesimal_sma, [p.sma])
        append!(planetesimal_Σ, [Σ])
    end
    protoplanet_sma = Number[]
    protoplanet_Σ = Number[]
    for i in 1:length(disk.protoplanets)
        if i < length(disk.protoplanets)
            outer_extent = disk.protoplanets[i+1].sma
        else
            outer_extent = 2disk.protoplanets[i].sma - disk.protoplanets[i-1].sma
        end
        Σ = mass2surfacedensity(disk.protoplanets[i].mass, disk.protoplanets[i].sma, outer_extent)
        append!(protoplanet_sma, [disk.protoplanets[i].sma])
        append!(protoplanet_Σ, [Σ])
    end

    if binned != false
        protoplanet_sma = [planetesimal_sma, disk.planetesimals[end].outer_extent]
        protoplanet_Σ = discrete_mass_to_surface_density(disk, protoplanet_sma)
    end

    mmsn = [10*sma^-1.5 for sma in planetesimal_sma]
    mmen = [100*sma^-1.5 for sma in planetesimal_sma]
    plot(layer(x=planetesimal_sma, y=planetesimal_Σ, Geom.line, color=["Planetesimals"]),
         layer(x=protoplanet_sma, y=protoplanet_Σ, Geom.point, color=["Protoplanets"]),
         layer(x=planetesimal_sma, y=mmsn, Geom.line, color=["MMSN"], Theme(default_color=color("red"))),
         layer(x=planetesimal_sma, y=mmen, Geom.line, color=["MMEN"], Theme(default_color=color("purple"))),
         Guide.xlabel("Semi-major axis [AU]"), Guide.ylabel("Surface density [g/cm^2]"), Guide.title(title),
         Scale.y_log10(minvalue=1e-2, maxvalue=1e5), Scale.x_log10,
         Theme(panel_fill=color("white"), panel_stroke=color("black"))
         )


    # df1 = DataFrame(x=planetesimal_sma, y=planetesimal_Σ, )


end


function discrete_mass_to_surface_density(disk, bin_bounds)
    sort!(disk)
    append!(bin_bounds, [1e100])
    current_bin = bin_bounds[1]
    current_bin_index = 1
    next_bin = bin_bounds[2]
    mass = zeros(length(bin_bounds)-1)
    for p in disk.protoplanets
        if p.sma >= current_bin
            while p.sma >= next_bin
                current_bin = next_bin
                current_bin_index += 1
                next_bin = bin_bounds[current_bin_index+1]
            end
            mass[current_bin_index] +=p.mass
        end
    end
    pop!(bin_bounds)
    pop!(mass)

    Σ = zeros(length(mass))
    for i=1:length(mass)
        Σ[i] = mass2surfacedensity(mass[i], bin_bounds[i], bin_bounds[i+1])
    end
    Σ
end




end






