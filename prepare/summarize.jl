"""
summarize: Energy and pseudologlikelihood summarization (prerequisite for CAMEOX)	
"""

using ArgParse, JLD, HDF5, GZip, StatsBase, LinearAlgebra

function mylogsumexp(x)
	X = maximum(x, dims = 2)
	return log.(sum(exp.(broadcast(-, x, X)), dims = 2)) + X
end

function convert_protein(prot_seq)
	prot_len = length(prot_seq)
	big_vector = zeros(Float32, prot_len * 21)
	#and fill it up like a very sparse row vector. One-hot.
	#Order is ARNDCQEGHILKMFPSTWYV-
	prot_fill = Dict{Char, Int64}('A'=>1, 'R'=>2, 'N'=>3, 'D'=>4, 'C'=>5, 'Q'=>6, 'E'=>7, 'G'=>8, 'H'=>9, 'I'=>10,
								  'L'=>11, 'K'=>12, 'M'=>13, 'F'=>14, 'P'=>15, 'S'=>16, 'T'=>17, 'W'=>18, 'Y'=>19, 'V'=>20,
                                  '-'=>21, 'X'=>21, '*'=>21, '.'=>21)  #added the dot, appearing in fasta from a3m msas
	chr_count = 0
	for chr in prot_seq
		big_vector[(21 * chr_count) + prot_fill[chr]] = Float32(1.0)
		chr_count += 1
	end
	return big_vector' #more convenient for it to be a row vector later.
end

function compute_energies(prot_mat, w1, w2)
	println("Computing energies...")
	flush(stdout)
	energies = prot_mat * w1 + diag(prot_mat * w2 * prot_mat')
	return energies #maybe this works?
end

function compute_psls(prot_mat, w1, w2, nNodes, nProt)
	println("Computing pseudolikelihoods...")
	flush(stdout)
	my_pv_w2 = prot_mat * w2
	my_ull = my_pv_w2 * prot_mat'
	my_pv_w1 = prot_mat * w1
	numerator = my_pv_w1 + diag(my_ull)
	partitions = reshape((w1' .+ my_pv_w2)', 21, nNodes * nProt)
	logZs = sum(reshape(mylogsumexp(partitions'), nNodes, nProt), dims = 1)'
	all_scores = numerator - logZs
	return -all_scores
end

function save_energies(computed_energies, gene_name, eng_path)
    local out_file
    try
	    out_file = open(joinpath(eng_path, "energy_$(gene_name).txt"), "w")
    catch err
        if isa(err, SystemError)
            out_file = open("energy_$(gene_name).txt", "w")
            println("WARNING! Saving MRF energy output to working directory")
        else
            rethrow()
        end
    finally
	    for energy in computed_energies
		    write(out_file, "$energy\n")
	    end
	    close(out_file)
    end
end

function save_psls(computed_psls, gene_name, pll_path)
    local out_file
    try
	    out_file = open(joinpath(pll_path, "psls_$(gene_name).txt"), "w")
    catch err
        if isa(err, SystemError)
            out_file = open("psls_$(gene_name).txt", "w")
            println("WARNING! Saving MRF pseudolikelihood output to working directory")
        else
            rethrow()
        end
    finally
	    for psl in computed_psls
		    write(out_file, "$psl\n")
	    end
	    close(out_file)
    end  
end

function run(gene_name, jld_file, msa_file, base_dir)
	println("Reading MSA file...")
	flush(stdout)
	in_file = open(msa_file)
	in_read = readlines(in_file)
	close(in_file)

	seqs = [uppercase(strip(xx)) for xx in filter(x -> x[1] != '>', in_read)]
	seqs = filter(seq -> !('B' in seq) && !('J' in seq) && !('O' in seq) && !('U' in seq) && !('Z' in seq), seqs) #try to avoid annoying sequences in alignment that have non-aas. X's are sorta okay.

	num_prot = length(seqs)
	num_aa = length(seqs[1]) #assume MSA is well formed.

	prot_mat = zeros(Float32, num_prot, num_aa * 21)
	i_count = 0
	for seq in seqs
		seq = uppercase(seq)
		i_count += 1
		prot_mat[i_count, :] = convert_protein(seq)
	end

	println("Reading JLD file...")
	flush(stdout)
	mrf = load(jld_file)
	w1 = mrf["w1"]
	w2 = mrf["w2"]

	@time computed_energies = compute_energies(prot_mat, w1, w2)
	@time computed_psls = compute_psls(prot_mat, w1, w2, num_aa, num_prot)

	#Write outputs to main directory.
	println("Computing done! Saving energies/psls...")
	flush(stdout)

    save_energies(computed_energies, gene_name, joinpath(base_dir, "energies/"))
    save_psls(computed_psls, gene_name, joinpath(base_dir, "psls/"))
end

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table s begin
        "prot_name"
		    help = "1st positional argument, name of gene"
            arg_type = String
        "prot_jld"
		    help = "2nd positional argument, path to JLD file"
            arg_type = String
		"prot_msa"
			help = "3rd positional argument, path to MSA (fasta) file"
            arg_type = String
		"base_dir"
			help = "4th positional argument, path to base dir of energy and psls directories"
            arg_type = String
	end
	return parse_args(s)
end

function summarize()
	println("=-= summarize = energy and pll =-= v0.2 - Jul 2023 =-= LLNL =-=")
	flush(stdout)
	
	# Parse arguments
	parsed_args = parse_commandline()
	prot_name = parsed_args["prot_name"]
	prot_jld = parsed_args["prot_jld"]
	prot_msa = parsed_args["prot_msa"]
	base_dir = parsed_args["base_dir"]

	if isfile(prot_msa)
		if isfile(prot_jld) && endswith(prot_jld, ".jld")
			run(prot_name, prot_jld, prot_msa, base_dir)
		else
			println("ERROR! $prot_jld is not a file or doesn't end with required .jld")
		end
	else
		println("ERROR! $prot_msa is not a file.")
	end
end

@time summarize()
