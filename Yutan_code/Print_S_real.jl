using DelimitedFiles,JLD

begin
    save_path = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/txt_results/"

    files = ["E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/Square_Spin_Ice_Measurement_tilted24sites_new_J2=1.0_hmin=0.01_hmax=1.5_Date_Mon_11_Oct_2021_22_42_19.jld"
    ]

    file = files[1]
    result = load(file, "result")
    hvals = result[1]

    outfile = save_path*"S_real_J2=$(result[2])"*".txt"

    open(outfile, "w") do io
        for (i, S_real) in enumerate(result[9][1])
            write(io, "h = $(hvals[i]) : \n")
            writedlm(io, S_real)
            write(io, "\n")
        end
    end

end
