export read_particle_data
export write_particles
export write_particles_minimal

"Reads the particle IDs and types of a given halo from the simulation data"
function read_halo_particles_IDs(haloID::String, simspecs::SimulationSpecs)::DataFrame
    input_particles_path = get_ahfparticles_filepath(simspecs)

    println("Reading particle IDs and types from $input_particles_path")

    halo_particle_IDs = Vector{Int64}()
    halo_particle_types = Vector{Int64}()

    open(input_particles_path, "r") do file
        readline(file) # Discard first line
        while !eof(file)
            line = readline(file)

            if endswith(line, haloID) # When the halo is found, extract N_particles and all the particle IDs
                println("Found specified halo")
                println("N_particles haloID")
                println(line)

                n_part = parse(Int, split(line)[1])

                sizehint!(halo_particle_IDs, n_part)
                sizehint!(halo_particle_types, n_part)

                @showprogress for i in 1:n_part
                    pID, ptype = parse.(Int, split(readline(file)))
                    push!(halo_particle_IDs, pID)
                    push!(halo_particle_types, ptype)
                end
                break # After finding the given halo and reading all the IDs, break the while loop
            end
        end
    end

    println()
    DataFrame(pid = halo_particle_IDs, ptype = halo_particle_types)
end
read_halo_particles_IDs(haloID::Int, simspecs::SimulationSpecs) = read_halo_particles_IDs(string(haloID), simspecs)

function read_particle_data_binary(halo_particles::DataFrame, simspecs::SimulationSpecs)
    particle_files = get_simparticle_filepaths(simspecs)
    pids_set = Set(halo_particles.pid) |> (ids -> convert.(UInt32, ids))
    println("Looking for $(size(halo_particles, 1)) DM particles")

    particles_dict = Dict(1 => DataFrame())

    read_blocksize(io::IO) = read(io, Int32)
    check_block_end(blocksize_init, blocksize_final) = blocksize_init == blocksize_final || error("Incorrect reading of data block.")

    for particle_file in particle_files
        if isempty(pids_set)
            println("All particles found, skipping the rest of the files")
            break
        end
        println("Reading $(split(particle_file, "/")[end])")

        # Binary file is read according the specification of GADGET-2 user guide section 6
        # https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf
        open(particle_file) do file
            file_dict = Dict()
            blocksize = read_blocksize(file)
            f_pos = position(file)
            # NOTE: Could make reading faster by only reading the high resolution particles
            # specified in n_particles[2] and then skipping the rest of the particles using
            # seek(file, f_pos + blocksize) and then proceeding to the next block
            n_particles = read!(file, Array{UInt32}(undef, 6))
            if n_particles[2] == 0
                println("No DM particles in this file")
                continue
            end
            seek(f, f_pos + blocksize)
            check_block_end(blocksize, read_blocksize(file))

            blocksize = read_blocksize(file)
            file_dict["Coordinates"] = read!(file, Array{Float32}(undef, 3, blocksize÷(4*3)))
            check_block_end(blocksize, read_blocksize(file))

            blocksize = read_blocksize(file)
            file_dict["Velocities"] = read!(file, Array{Float32}(undef, 3, blocksize÷(4*3)))
            check_block_end(blocksize, read_blocksize(file))

            blocksize = read_blocksize(file)
            file_dict["ParticleIDs"] = read!(file, Array{UInt32}(undef, blocksize÷4))
            check_block_end(blocksize, read_blocksize(file))

            blocksize = read_blocksize(file)
            file_dict["Masses"] = read!(file, Array{Float32}(undef, blocksize÷4))
            check_block_end(blocksize, read_blocksize(file))
        end

        mask = zeros(Bool, length(file_dict["ParticleIDs"]))
        tmap!(in(pids_set), mask, file_dict["ParticleIDs"])
        println("Found $(sum(mask)) DM particles")
        setdiff!(pids_set, file_dict["ParticleIDs"][mask]) # Remove the found particles

        update_particles_dict!(particles_dict, file_dict, mask, 1)
    end

    println("Extracted $(size(particles_dict[1], 1)) DM particles")
    println("Difference from expected: $(size(particles_dict[1], 1) - size(halo_particles, 1))")

    println()
    convert_particles_dict_to_df(particles_dict)
end

function read_particle_data_hdf5(halo_particles::DataFrame, simspecs::SimulationSpecs)
    particle_types = [0, 1, 4, 5]
    particles_dict = Dict(type => DataFrame() for type in particle_types)
    halo_particles_grouped = groupby(halo_particles, :ptype)

    println("In halo particleIDs file:")
    for particle_type in particle_types
        pcount = count(==(particle_type), halo_particles.ptype)
        println("Type $(particle_type) has $(pcount) particles")
        if pcount == 0
            deleteat!(particle_types, findfirst(==(particle_type), particle_types))
        end
    end
    println()

    particle_files = get_simparticle_filepaths(simspecs)

    for particle_file in particle_files
        particles = h5open(particle_file, "r")
        println("Reading $(split(particle_file, "/")[end])")
        println("Matches:")

        for particle_type in particle_types
            halo_ids = halo_particles_grouped[(ptype = particle_type,)].pid |> (ids -> convert.(UInt32, ids)) |> Set
            try # Try and catch if the file doesn't have a certain particle type (like black holes)
                particles["PartType$(particle_type)"]
            catch e
                if isa(e, KeyError)
                    println("Type $(particle_type): 0 particles")
                    continue
                else
                    rethrow(e)
                end
            end

            all_ids = read(particles, "PartType$(particle_type)/ParticleIDs")
            mask = zeros(Bool, length(all_ids))
            tmap!(in(halo_ids), mask, all_ids)
            println("Type $(particle_type): $(sum(mask)) particles")

            group_dict = read(particles, "PartType" * string(particle_type))
            update_particles_dict!(particles_dict, group_dict, mask, particle_type)
        end

        close(particles)
    end

    println()
    println("Extracted:")
    for particle_type in particle_types
        print("Type $(particle_type) has $(size(particles_dict[particle_type], 1)) particles ")
        println("(diff $(size(particles_dict[particle_type], 1) - count(==(particle_type), halo_particles.ptype)))")
    end

    println()
    convert_particles_dict_to_df(particles_dict)
end
read_particle_data(halo_particles::DataFrame, simspecs::SimulationSpecs) = simspecs.galfor ? read_particle_data_hdf5(halo_particles, simspecs) : read_particle_data_binary(halo_particles, simspecs)
read_particle_data(haloID::Union{Int, String}, simspecs::SimulationSpecs) = read_particle_data(read_halo_particles_IDs(haloID, simspecs), simspecs)

function update_particles_dict!(particles_dict::Dict{Int64, DataFrame}, group_dict::Dict, mask::Vector{Bool}, particle_type::Int64)
    # Separate the matrix datasets into individual columns with numbers appended to the property key
    for prop_key in keys(group_dict)
        if typeof(group_dict[prop_key]) <: Matrix
            for i in 1:size(group_dict[prop_key], 1)
                group_dict[prop_key * string(i)] = group_dict[prop_key][i, :]
            end
            delete!(group_dict, prop_key)
        end
    end

    particles_dict[particle_type] = vcat(
        particles_dict[particle_type],
        DataFrame(group_dict)[mask, :] # Only write the lines that are in the mask
    )
end

function convert_particles_dict_to_df(particles_dict::Dict{Int64, DataFrame})
    particles_df = DataFrame()
    particle_types_dict = Dict(
        0 => :gas,
        1 => :dm,
        4 => :stars,
        5 => :bh,
    )

    for ptype_num in keys(particles_dict)
        temp_df = particles_dict[ptype_num]
        length_df = size(temp_df, 1)
        temp_df.ptype = fill(particle_types_dict[ptype_num], length_df)
        append!(particles_df, temp_df; cols = :union)
    end

    if !(:GFM_StellarFormationTime in names(particles_df))
        particles_df.GFM_StellarFormationTime = fill(missing, size(particles_df, 1))
    end

    particles_df
end

function write_particles(particles_df::DataFrame, haloID::Union{Int, String}, simspecs::SimulationSpecs)
    fout = "HESTIA_$(simspecs.simID)_$(simspecs.n_particles)_halo$(haloID)_particles.arrow"
    Arrow.write(joinpath(data_dir, "haloparticles", fout), particles_df) # data_dir is a const in utils.jl
end
write_particles(haloID::Union{Int, String}, simspecs::SimulationSpecs) = write_particles(read_particle_data(haloID, simspecs), haloID, simspecs)

function write_particles_minimal(particles_df::DataFrame, haloID::Union{Int, String}, simspecs::SimulationSpecs)
    particles_df_minimal = select(particles_df,
        :ParticleIDs,
        :ptype,
        :Coordinates1,
        :Coordinates2,
        :Coordinates3,
        :Velocities1,
        :Velocities2,
        :Velocities3,
        :Masses,
        :GFM_StellarFormationTime,
    )
    fout = "HESTIA_$(simspecs.simID)_$(simspecs.n_particles)_halo$(haloID)_minimal_particles.arrow"
    Arrow.write(joinpath(data_dir, "haloparticles", fout), particles_df_minimal) # data_dir is a const in utils.jl
end
write_particles_minimal(haloID::Union{Int, String}, simspecs::SimulationSpecs) = write_particles_minimal(read_particle_data(haloID, simspecs), haloID, simspecs)

