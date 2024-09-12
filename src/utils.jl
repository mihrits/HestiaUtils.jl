
export SimulationSpecs

include("snapshot2z.jl") # Dictionary that converts snapshot::Int to correct snapshot-redshift filename string

const base_dir = "/store/clues/HESTIA"
const project_dir = "~/projects/dmannih"

"Struct containing the simulation specifications"
struct SimulationSpecs
    simID::String
    snapshot::Int64
    n_particles::Int64

    function SimulationSpecs(simID::String, snapshot::Int64 = 127, n_particles::Int64 = 8192)
        new(simID, snapshot, n_particles)
    end
end

"Returns the base filepath of the AHF files for a given simulation"
function get_ahfbasepath(simspecs::SimulationSpecs)::String
    # 8192 particle simulations AHF outputs include "2x2.5Mpc" in the directory name
    AHF_output_dir = simspecs.n_particles == 8192 ? "AHF_output_2x2.5Mpc" : "AHF_output"

    joinpath(
        base_dir, # const in utils.jl
        "RE_SIMS",
        $(simsspecs.n_particles),
        "GAL_FOR",
        $(simspecs.simID),
        AHF_output_dir,
        "HESTIA_100Mpc_$(simspecs.n_particles)_$(simspecs.simID).$(snapshot2z[simspecs.snapshot]).AHF_",
    )
end

"Returns the filepath of the AHF particles file for a given simulation"
function get_ahfparticles_filepath(simspecs::SimulationSpecs)::String
    get_ahfbasepath(simspecs, snapshot2z, base_dir) * "particles"
end

"Returns the filepath of the AHF profiles file for a given simulation"
function get_ahfprofiles_filepath(simspecs::SimulationSpecs)::String
    get_ahfbasepath(simspecs, snapshot2z, base_dir) * "profiles"
end

"Returns the filepath of the AHF halos file for a given simulation"
function get_ahfhalos_filepath(simspecs::SimulationSpecs)::String
    get_ahfbasepath(simspecs, snapshot2z, base_dir) * "halos"
end

function get_simparticle_filepaths(simspecs::SimulationSpecs)::Vector{String}
    # 8192 particle simulations outputs include "2x2.5Mpc" in the directory name
    output_dir = simspecs.n_particles == 8192 ? "output_2x2.5Mpc" : "output"

    snapdir_path = joinpath(
        base_dir, # const in utils.jl
        "RE_SIMS",
        $(simspecs.n_particles),
        "GAL_FOR",
        $(simspecs.simID),
        output_dir,
        "snapdir_$(simspecs.snapshot)",
    )

    readdir(snapdir_path, join = true, sort = true)
end

