"""
raw2jld: Convert CCMpred RAW file to JLD file (prerequisite for CAMEOX)	
"""

using ArgParse, JLD, HDF5

function sub_read(w2::Array{Float64, 2}, co1::Int64, co2::Int64, cur_rows::Array{Float64, 2})
	w2[1 + co1 * 21: (co1 + 1) * 21, 1 + co2 * 21: (co2 + 1) * 21] = cur_rows[1:end, 1:end]
	w2[1 + co2 * 21: (co2 + 1) * 21, 1 + co1 * 21: (co1 + 1) * 21] = cur_rows[1:end, 1:end]'
end

function read_raw(raw_file, out_file)
	in_file = open(raw_file)
	raw_read = readlines(in_file)
	close(in_file)

	ncol = findfirst(n -> n[1] == '#', raw_read) - 1

	w1 = zeros(ncol, 20)
	count = 0
	for line in raw_read[1:ncol]
		count += 1
		w1[count, :] = map(xx -> parse(Float64, xx), split(strip(line), "\t"))
	end

	w2 = zeros(ncol * 21, ncol * 21)
	co1 = 0
	co2 = 0
	row = 0
	cur_rows = zeros(21, 21)
	for line in raw_read[(ncol + 1):(end-1)]
		if line[1] == '#'
			if !(co1 == 0 && co2 == 0)
				sub_read(w2, co1, co2, cur_rows)
			end

			row = 0
			pos = split(strip(line))
			cur_rows *= 0.0 #reset
			co1 = parse(Int64, pos[2])
			co2 = parse(Int64, pos[3])
		else
			row += 1
			cur_rows[row, :] = map(mm -> parse(Float64, mm), split(strip(line), "\t"))
		end
	end

	w1_x = hcat(w1, zeros(ncol, 1))'
	w1_alt = reshape(w1_x, (21 * ncol, 1))

	save(out_file, "w1", Float64.(w1_alt), "w2", w2)
	println("raw2jld OK! Saved $out_file")

end

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table s begin
		"raw_file"
			help = "1st positional argument, path to CCMpred RAW file"
            arg_type = String
			required = true
        "jld_file"
		    help = "2nd positional argument, path to JLD file"
            arg_type = String
			required = true
	end
	return parse_args(s)
end

function raw2jld()
	println("=-= raw2jld = CCMpred RAW to JLD =-= v0.2 - Aug 2023 =-= LLNL =-=")
	flush(stdout)
	
	# Parse arguments
	parsed_args = parse_commandline()
	raw_file = parsed_args["raw_file"]
	jld_file = parsed_args["jld_file"]

	# Checks and main call
	if isfile(raw_file)
		if endswith(raw_file, ".raw")
			if endswith(jld_file, ".jld")
				println("raw2jld converting $raw_file... Please wait...")
				read_raw(raw_file, jld_file)
			else
				println("ERROR! Second file should end in .jld")
			end
		else
			println("ERROR! First file should end in .raw")
		end
	else
		println("ERROR! Incorrect input, check if raw file exists.")
	end
end

@time raw2jld()
