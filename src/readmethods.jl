export read_ahfmergertree

function read_ahfmergertree(haloID::Int, simspecs::SimulationSpecs)
    filein = get_ahfmergertree_filepath(haloID, simspecs)
    if !isfile(filein) return error("File $(filein) does not exist.") end

    # The first column is delimited by " " and all the others by "\t"
	# Only read the first 4 columns (z and ID are read as one column, because of wrong delim)
    mergertree_readin_opts = (comment = "#", header = false, delim = "\t", ignorerepeated = false, select = 1:12)

    mergertree_header = (
        "Column2"  => :haloID,
        "Column4"  => :Mvir,
        "Column6"  => :center_x,
        "Column7"  => :center_y,
        "Column8"  => :center_z,
        "Column9"  => :center_vx,
        "Column10" => :center_vy,
        "Column11" => :center_vz,
        "Column12" => :Rvir,
    )
    @chain begin
        CSV.read(filein, DataFrame; mergertree_readin_opts...)
        # Separate z and ID and parse them into Float and Int
        select("Column1" => ByRow(split) => [:z, :ID], mergertree_header...)
        transform!(:z => ByRow(z -> parse(Float64, z)) => :z, :ID => ByRow(id -> parse(Int, id)) => :ID)
    end
end