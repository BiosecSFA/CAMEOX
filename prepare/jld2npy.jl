"""
jld2npy: convert CAMEOX JLD to numpy format
"""
# Based on CAMEOS' convert_jld_to_npy.jl

using NPZ, JLD

function jld2npy()
    println("=-= CAMEOX jld2npy.jl =-= v0.1 - Ago 2021 =-= by LLNL =-=")    
	in_file = ARGS[end]
	if isfile(in_file)
		if endswith(in_file, ".jld")
			model = load(in_file)

			w1 = model["w1"]
			w2 = model["w2"]

			npzwrite(replace(in_file, ".jld" => "_w1") * ".npy", w1)
			npzwrite(replace(in_file, ".jld" => "_w2") * ".npy", w2)
		else
			error("$in_file needs to be a JLD file and have suffix .jld")
		end
	else
		error("$in_file is not a file!")
	end
    filenames = replace(in_file, ".jld" => "_w[1-2]") * ".npy"
    println("OK! Generated $filenames files")
end

@time jld2npy()
