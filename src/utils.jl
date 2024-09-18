
export SimulationSpecs
export get_galvectors_from_profile

include("snapshot2z.jl") # Dictionary that converts snapshot::Int to correct snapshot-redshift filename string

const hestia_dir = "/store/clues/HESTIA"
const data_dir = "/store/erebos/moortis/projects/dmannih/data"

"Struct containing the simulation specifications"
struct SimulationSpecs
    simID::String
    snapshot::Int64
    n_particles::Int64

    function SimulationSpecs(simID::String, snapshot::Int64 = 127, n_particles::Int64 = 8192)
        new(simID, snapshot, n_particles)
    end
end

function check_haloID_simspecs_compatibility(haloID::Int, simspecs::SimulationSpecs)::Nothing
    if evalpoly(10, digits(haloID)[13:end]) != simspecs.snapshot
        error("HaloID $(haloID) is not compatible with simulation snapshot $(simspecs.snapshot).")
    end
end

"Returns the base filepath of the AHF files for a given simulation"
function get_ahfbasepath(simspecs::SimulationSpecs)::String
    # 8192 particle simulations AHF outputs include "2x2.5Mpc" in the directory name
    AHF_output_dir = simspecs.n_particles == 8192 ? "AHF_output_2x2.5Mpc" : "AHF_output"

    joinpath(
        hestia_dir, # const in utils.jl
        "RE_SIMS",
        string(simspecs.n_particles),
        "GAL_FOR",
        simspecs.simID,
        AHF_output_dir,
        "HESTIA_100Mpc_$(simspecs.n_particles)_$(simspecs.simID).$(snapshot2z_dict[simspecs.snapshot]).AHF_",
    )
end

"Returns the filepath of the AHF particles file for a given simulation"
function get_ahfparticles_filepath(simspecs::SimulationSpecs)::String
    get_ahfbasepath(simspecs) * "particles"
end

"Returns the filepath of the AHF profiles file for a given simulation"
function get_ahfprofiles_filepath(simspecs::SimulationSpecs)::String
    get_ahfbasepath(simspecs) * "profiles"
end

"Returns the filepath of the AHF halos file for a given simulation"
function get_ahfhalos_filepath(simspecs::SimulationSpecs)::String
    get_ahfbasepath(simspecs) * "halos"
end

function get_simparticle_filepaths(simspecs::SimulationSpecs)::Vector{String}
    # 8192 particle simulations outputs include "2x2.5Mpc" in the directory name
    output_dir = simspecs.n_particles == 8192 ? "output_2x2.5Mpc" : "output"

    snapdir_path = joinpath(
        hestia_dir, # const in utils.jl
        "RE_SIMS",
        string(simspecs.n_particles),
        "GAL_FOR",
        simspecs.simID,
        output_dir,
        "snapdir_$(lpad(simspecs.snapshot, 3, '0'))",
)

    filter!(endswith(".hdf5"), readdir(snapdir_path, join = true, sort = true))
end

function get_ahfmergertree_filepath(haloID::Int, simspecs::SimulationSpecs)::String
    check_haloID_simspecs_compatibility(haloID, simspecs)
    # 8192 particle simulations AHF outputs include "2x2.5Mpc" in the directory name
    AHF_output_dir = simspecs.n_particles == 8192 ? "AHF_output_2x2.5Mpc" : "AHF_output"

    joinpath(
        hestia_dir, # const in utils.jl
        "RE_SIMS",
        string(simspecs.n_particles),
        "GAL_FOR",
        simspecs.simID,
        AHF_output_dir,
        "HESTIA_100Mpc_$(simspecs.n_particles)_$(simspecs.simID).$(simspecs.snapshot)_halo_$(haloID).dat",
    )
end

"""
    get_galvectors_from_profile(profile::DataFrame, r_min::Real)

Returns the Ec (normal of the disk) and Ea (major disk axis) vectors from a given AHF
profile DataFrame at the smallest radius that is r >= r_min [kpc].
"""
function get_galvectors_from_profile(profile::DataFrame, r_min::Real)
    idx = findfirst(>=(r_min), profile.r)
    collect(profile[idx, [:Ecx, :Ecy, :Ecz]]), collect(profile[idx, [:Eax, :Eay, :Eaz]])
end

"""
    get_galvectors_from_profile(haloID::Int, simspec::SimulationSpecs, r_min::Real)

Returns the Ec (normal of the disk) and Ea (major disk axis) vectors for `haloID` and
simulation `simspecs` at the smallest radius that is r >= r_min [kpc].
"""
function get_galvectors_from_profile(haloID::Int, simspec::SimulationSpecs, r_min::Real)
    get_galvectors_from_profile(read_ahfprofile(haloID, simspec), r_min)
end