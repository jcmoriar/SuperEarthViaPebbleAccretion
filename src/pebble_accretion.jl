using MKS
using Roots


############################################################################
# General Functions
############################################################################

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


function gasdensity(a, t, β_0=500)
    """Calculate gas density as a function of semi-major axis.

        Args:
            a: Semi-major axis (AU).
            t: Time in years.
        Return:
            Gas density in g/cm^3.

    """

    scale_height = 0.033(a^1.25) * AU * a * 100  # cm
    sigma = β_0 * exp(-t/3e6) / a  # g/cm^2
    sigma / scale_height / 2
end


function gassurfacedensity(a, t, β_0=500)
    """Calculate gas surface density as a function of semi-major axis.

        Args:
            a: Semi-major axis (AU).
            t: Time in years.
        Return:
            Gas density in g/cm^3.

    """

    β_0 * exp(-t/3e6) / a  # g/cm^2
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


function surfacedensity2mass(Σ, inner, outer)
    """ Calculate the mass in an annulus.

        Args:
            Σ: Surface density in annulus [g/cm^2].
            inner: Inner edge of annulus in AU.
            outer: Outer edge of annulus in AU.

        Return:
            Mass in annulus in M_Earth.

    """
    area = π*((outer*100AU)^2 - (inner*100AU)^2)
    area * Σ / 1000MEARTH
end


############################################################################
# Pebble Production and Surface Density
############################################################################

function pebbleproduction(t, β=500., m_star=1., Z=0.01)
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
    9.5e-5(β*exp(-t/3e6)/500.) * m_star^(1./3.) * (Z/0.01)^(5./3.) * (t/1e6)^(-1./3.)
end


function pebbledensity(a, t; r0=1, ρ=2, pebble_flux=nothing)
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
    τ_dis = 3e6
    β = 500 * exp(-t/τ_dis)
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
    ρ = 2
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

function alpha_p(a, ρ_s=3):
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


function accretion_regime(a, t, r_p, r_0; ρ_pl=3, ρ_peb=2)
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


function impact_parameter(a, t, r_p, r_0; ρ_pl=3, ρ_peb=2)
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

function relative_velocity(a, b, r0, t; α_t=0.01, e=0)
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


function growth_rate(a, t, r_p, r0; e=0, α_t=0.0001, ρ_pl=3, ρ_peb=2, pebble_flux=nothing)
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
        individual_mass = sphere_mass(body_size, 3.) / MEARTH
        num_bodies = total_mass/individual_mass
        # println(num_bodies)
    else
        individual_mass = total_mass/num_bodies
        body_size = sphere_radius(individual_mass, 3.)
    end
    Ṁ_individual = growth_rate(a, t, body_size, 1, pebble_flux=pebble_flux)
    # println(num_bodies, " ", individual_mass/Ṁ_individual, " ", pebble_flux)
    Ṁ_individual * num_bodies
end


#######################################################################
# Full Simulation
#######################################################################

type Disk
    annuli::Array
    annuli_bounds::Array
    surface_density::Array
    max_surface_density::FloatingPoint
    embryos::Array
end

type Embryo
    sma::FloatingPoint
    mass::FloatingPoint
    eccentricity::FloatingPoint
end


function Disk(inner_edge::FloatingPoint, outer_edge::FloatingPoint, num_annuli::Integer, max_surface_density)
    """ Constructor for type PlanetesimalDisk.

        Args:
            inner_edge: The inner edge of the plnaetesimal disk in AU.
            outer_edge: The outer edge of the plnaetesimal disk in AU.
            num_annuli: Defines resolution of simulation.
            max_surface_density: The threshold at which additional mass
                is converted to a new planetary embryo.

        Returns:
            A new PlanetesimalDisk object.

    """
    annuli = [10^loga for loga in linspace(log10(inner_edge), log10(outer_edge), num_annuli+1)]
    embryos = Embryo[]
    annuli_bounds = []
    for i = 1:num_annuli
        annuli_bounds = [annuli_bounds, annuli[i], annuli[i+1]]
        append!(embryos, [Embryo(annuli[i], surfacedensity2mass(3., annuli[i], annuli[i+1]), 0)])
    end

    annuli_bounds = transpose(reshape(annuli_bounds, 2, num_annuli))


    ########
    # I'm not sure what to do about the initial surface density.
    # I'll start with a constant 3 g/cm^3 and go from there
    ########
    surface_density = zeros(num_annuli)+3
    # embryos = Embryo[]
    # for i=1:100
    #     append!(embryos, [Embryo(rand()*(outer_edge-inner_edge)+inner_edge,1e-2,0)])
    # end
    Disk(annuli, annuli_bounds, surface_density, max_surface_density, embryos) #random embryos for now
end


function Disk(inner_edge::FloatingPoint, outer_edge::FloatingPoint,
              total_mass::FloatingPoint, min_sep::FloatingPoint,
              max_sep::FloatingPoint)
    """ Constructs an EmbryoDisk with the provided parameters.

        Args:
            inner_edge: Inner extent of disk in AU.
            outer_edge: Outer extent of disk in AU.
            total_mass: Total mass of embryos in Earth masses.
            min_sep: Minimum separation between embryos in mutual hill radii.
            max_sep: Maximum separation between embryos in mutual hill radii.

    """
    x = "I'll leave this one for later"
end


function accrete_pebbles!(disk::Disk, time, dt)
    """ Calculate the pebble accretion rate at each annulus and apply it for the provided time interval.

        Args:
            disk: <Disk> The protoplanetary disk for which to calculate the pebble accretion rates.
            time: The time in years.
            dt: The time step in years.

    """
    # Start by determining which annulus each embryo lies in
    sortkey(embryo) = embryo.sma
    sort!(disk.embryos, by=sortkey)
    embryo_annuli=Dict()
    current_annulus = 0
    next_annulus = disk.annuli[1]
    next_index = 1
    max_index = length(disk.annuli)
    embryo_annuli[0] = Embryo[]
    plotquantity = Number[]
    quantityname = "Embryo surface density [g/cm^2]"
    for embryo in [disk.embryos, Embryo(1e50, 1, 0)] # one extra distant embryo to ensure all annuli are initialized
        while embryo.sma >= next_annulus
            current_annulus=next_annulus
            next_index += 1
            embryo_annuli[current_annulus] = Embryo[]
            if next_index > max_index
                next_annulus=1e100
            else
                next_annulus = disk.annuli[next_index]
            end
        end
        append!(embryo_annuli[current_annulus], [embryo])
    end
    pop!(embryo_annuli[disk.annuli[end]])  # Remove the extra distant embryo
    # for key in sort(collect(keys(embryo_annuli)))
    #     println(key, embryo_annuli[key])
    # end

    # Now determine and apply the pebble accretion rate to bodies in each annulus
    pebble_flux = pebbleproduction(time)
    pebble_flux -= accrete_pebbles!(embryo_annuli[disk.annuli[end]], time, dt, pebble_flux)
    total_accreted = 0
    for i=length(disk.surface_density):-1:1
        # println(pebble_flux)
        embryo_sink = accrete_pebbles!(embryo_annuli[disk.annuli[i]], time, dt, pebble_flux)
        massinannulus = surfacedensity2mass(disk.surface_density[i], disk.annuli_bounds[i,:]...)
        # println(disk.annuli_bounds[i,1], " ", disk.annuli_bounds[i,2], " ", massinannulus)
        # println(pebble_flux)
        planetesimal_sink = total_accretion_rate(disk.annuli[i], time, massinannulus,
                                                 body_size=100., pebble_flux=pebble_flux)
        # println(disk.annuli_bounds[i,1], " ", disk.annuli_bounds[i,2], " ", planetesimal_sink)
        embryomass=0
        for e in embryo_annuli[disk.annuli[i]]
            embryomass += e.mass
        end
        append!(plotquantity, [mass2surfacedensity(embryomass+1e-20, disk.annuli_bounds[i,:]...)])
        Δmass = planetesimal_sink * dt
        ΔΣ_planetesimal = mass2surfacedensity(Δmass, disk.annuli_bounds[i,:]...)
        disk.surface_density[i] += ΔΣ_planetesimal
        if (embryo_sink + planetesimal_sink) > pebble_flux
            println("The accretion rate is higher than the flux at ", disk.annuli[i], " AU")
            # println("Warning: More than 5% of pebbles are being accreted in a single annulus. This may affect the calculation. ", disk.annuli_bounds[i,1])
            # println((embryo_sink + planetesimal_sink)/pebble_flux)
        end
        pebble_flux -= embryo_sink + planetesimal_sink
        total_accreted += (embryo_sink + planetesimal_sink) * dt
    end
    convert_planetesimals_to_embryos!(disk)
    # println("Total Accreted mass: ", total_accreted/pebbleproduction(time)/dt)
    (reverse!(plotquantity), quantityname)
end


function accrete_pebbles!(embryos::Array{Embryo}, time, dt, pebble_flux)
    """ Calculate and apply the pebble accretion rate to an array of embryos.

    Args:
        embryos: Array of embryos.
        time: Time for calculation in years.
        dt: Timestep in years.
        pebble_flux: The flux of pebbles available to accrete.

    Return:
        The combined accretion rate of the embryos in M_Earth/yr.

    """
    total_accretion=0
    for embryo in embryos
        Σ_peb = pebbledensity(embryo.sma, time, pebble_flux=pebble_flux)
        radius = sphere_radius(embryo.mass, 3.)
        Ṁ = growth_rate(embryo.sma, time, radius,
                        1., pebble_flux=pebble_flux) # 1cm particles; 3g/cm^3 planet density
        # if Ṁ < 0
        #     println(embryo.sma)
        #     println(radius)
        #     println(time)
        #     println(pebble_flux)
        # end
        embryo.mass += Ṁ * dt
        total_accretion += Ṁ
    end
    total_accretion
end


function convert_planetesimals_to_embryos!(disk::Disk)
    # If there is too much mass in planetesimals create a new embryo
    accumulated_mass = 0
    inner_bound = disk.annuli_bounds[1,1]
    included_annuli = Integer[]
    excess_Σ = FloatingPoint[]
    for i in 1:length(disk.annuli)-1
        if disk.surface_density[i] > disk.max_surface_density
            append!(included_annuli, [i])
            append!(excess_Σ, [disk.surface_density[i]-disk.max_surface_density])
            accumulated_mass += surfacedensity2mass(disk.surface_density[i]-disk.max_surface_density, disk.annuli_bounds[i,:]...)
            if accumulated_mass > 1e-2
                num_bodies = iceil(accumulated_mass/1e-1)
                for j in 1:num_bodies
                    new_embryo = Embryo(rand()*(disk.annuli_bounds[i,2] - inner_bound)+inner_bound, accumulated_mass/num_bodies, 0)
                    append!(disk.embryos, [new_embryo])
                end
                disk.surface_density[included_annuli] -= excess_Σ
                included_annuli = Integer[]
                excess_Σ = FloatingPoint[]
                inner_bound = disk.annuli_bounds[i,2]
            end
        end
    end
end













