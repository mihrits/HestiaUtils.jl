export read_ahfmergertree
export read_ahfprofile
export read_ahfhalos
export read_ahfsubhalos

function read_ahfmergertree(haloID::Int, simspecs::SimulationSpecs)
    check_haloID_simspecs_compatibility(haloID, simspecs)

    filein = get_ahfmergertree_filepath(haloID, simspecs)
    if !isfile(filein) return error("File $(filein) does not exist.") end

    # The first column is delimited by " " and all the others by "\t"
	# Only read the first 4 columns (z and ID are read as one column, because of wrong delim)
    mergertree_readin_opts = (comment = "#", header = false, delim = "\t", ignorerepeated = false, select = 1:12)

    mergertree_header = (
        "Column2"  => :hostID,
        "Column4"  => :Mvir,
        "Column6"  => :Xc,
        "Column7"  => :Yc,
        "Column8"  => :Zc,
        "Column9"  => :VXc,
        "Column10" => :VYc,
        "Column11" => :VZc,
        "Column12" => :Rvir,
    )
    @chain begin
        CSV.read(filein, DataFrame; mergertree_readin_opts...)
        # Separate z and ID and parse them into Float and Int
        select("Column1" => ByRow(split) => [:z, :haloID], mergertree_header...)
        transform!(:z => ByRow(z -> parse(Float64, z)) => :z, :haloID => ByRow(id -> parse(Int, id)) => :haloID)
    end
end

function read_ahfprofile(haloID::Int, simspecs::SimulationSpecs)
    check_haloID_simspecs_compatibility(haloID, simspecs)

    filein = get_ahfprofiles_filepath(simspecs)
    if !isfile(filein) return error("File $(filein) does not exist.") end

    target_halo_nr = parse(Int, string(haloID)[4:end])
    halo_profile_lines = Vector{String}()
    header = Vector{Symbol}()

    open(filein, "r") do file
        header = [Symbol(h) for h in split(readline(file)[2:end])]

        current_halo_nr = 1
        line = readline(file)
        r1 = abs(parse(Float64, split(line)[1]))

        if current_halo_nr == target_halo_nr
            push!(halo_profile_lines, line)
        end

        while !eof(file)
            line = readline(file)
            r2 = abs(parse(Float64, split(line)[1]))
            if r2 < r1
                current_halo_nr += 1
            end
            if current_halo_nr == target_halo_nr
                push!(halo_profile_lines, line)
            elseif current_halo_nr > target_halo_nr
                break
            end
            r1 = r2
        end
    end
    halo_profile = CSV.read(
        IOBuffer(join(halo_profile_lines, "\n")),
        DataFrame; header = header
    )
    select!(halo_profile,
        "r(1)" => :r,
        "npart(2)" => :npart,
        "M_in_r(3)" => :M_in_r,
        "dens(5)" => :dens,
        "Eax(14)" => :Eax,
        "Eay(15)" => :Eay,
        "Eaz(16)" => :Eaz,
        "Ebx(17)" => :Ebx,
        "Eby(18)" => :Eby,
        "Ebz(19)" => :Ebz,
        "Ecx(20)" => :Ecx,
        "Ecy(21)" => :Ecy,
        "Ecz(22)" => :Ecz,
        "b(12)" => :b,
        "c(13)" => :c,
    )
end

function read_ahfhalos(simspecs::SimulationSpecs)
    filein = get_ahfhalos_filepath(simspecs)
    if !isfile(filein) return error("File $(filein) does not exist.") end

    halos_readin_opts = (comment = "#", header = false, delim = "\t", ignorerepeated = false, select = 1:12)

    halos_header = (
        "Column1"  => :haloID,
        "Column2"  => :hostID,
        "Column3"  => :numSubStruct,
        "Column4"  => :Mvir,
        "Column6"  => :Xc,
        "Column7"  => :Yc,
        "Column8"  => :Zc,
        "Column9"  => :VXc,
        "Column10" => :VYc,
        "Column11" => :VZc,
        "Column12" => :Rvir,
    )

    select!(CSV.read(filein, DataFrame; halos_readin_opts...), halos_header...)
end

function read_ahfsubhalos(haloID::Int, path_mtree::String; verbose = true)
    if !isfile(path_mtree)
        return error("File $(path_mtree) does not exist.")
    end

    subhaloIDs = Vector{Int64}()

    open(path_mtree, "r") do file
        readline(file) # Skip initial line with nr of host halos
        while !eof(file)
            ID, N_subhalos = readline(file) |> split |> x -> parse.(Int, x)

            if ID == haloID
                if verbose
                    println("Found specified halo with ID $haloID")
                    println("N_subhalos: $N_subhalos")
                end

                sizehint!(subhaloIDs, N_subhalos-1)
                readline(file) # Skip the first ID because it is the haloID itself

                for _ in 1:N_subhalos-1
                    subhaloID = parse(Int, readline(file))
                    push!(subhaloIDs, subhaloID)
                end

                break
            else
                for _ in N_subhalos
                    readline(file) # Skip the subhalo IDs of this halo
                end
            end
        end
    end

    if verbose println() end
    subhaloIDs
end

function read_ahfsubhalos(haloID::Int, simspecs::SimulationSpecs; verbose = true)
    path_mtree = get_ahfsubhalos_filepath(simspecs)
    read_ahfsubhalos(haloID, path_mtree; verbose = verbose)
end