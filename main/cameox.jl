"""
CAMEOX : CAMEOs eXtended main module
"""

module CAMEOX

include("utils.jl")
include("types.jl")
include("lookup.jl")
include("std_setup.jl")
include("math.jl")
include("optimize.jl")
include("mrf.jl")
include("bio_seq.jl")

using ArgParse, Dates, Distributions, FileIO, JLD, Logging
using Printf, Random, Statistics, StatsBase
#using Gtk, Profile, ProfileView

"""
Aux: Get a JLD saveable population of variants
"""
function population_to_saveable(cur_pop;
                                deg_wt_apll::Float32=0.0f0,
								mark_wt_apll::Float32=0.0f0)
	new_pop = types.SaveChrome[]

	for indiv in cur_pop
		push!(
            new_pop,
            types.SaveChrome(
                indiv.path, indiv.full_sequence, indiv.deg_nuc, indiv.deg_map,
                indiv.deg_trns, indiv.deg_trne, indiv.deg_d, indiv.deg_skip,
                indiv.deg_insert, indiv.mark_nuc, indiv.mark_map,
                indiv.mark_trns, indiv.mark_trne, indiv.mark_d,
                indiv.mark_skip, indiv.mark_insert,	indiv.deg_prob,
                indiv.deg_base_E, indiv.deg_seq, indiv.deg_prob-deg_wt_apll, 
				indiv.mark_prob, indiv.mark_base_E, indiv.mark_seq,
                indiv.mark_prob-mark_wt_apll, indiv.first_weight))
	end
	return new_pop
end

"""
set_up_and_optimize is the code that sets up appropriate information about genes we want to target for double-encoding and
orchestrates the initial HMM double-encoding solution and then iteratively improves on that with greedy pseudolikelihood improvements.
Code also logs information through Logging package.
Genes are often referred to as mark/deg. This is historic from times of thinking of "marker" genes and "designated essential gene".
At this point this distinction is not meaningful but is kept to avoid equally arbitrary x_name, y_name conventions.
Arguments:
	log_io: IO object for logging 
	rand_barcode: barcode for this entanglement pair
	paths: Paths object with relevant paths
	mark_name : gene ID for 'mark' gene. Needed as a key for looking up some values associated with genes in text files.
	deg_name : gene ID for 'deg' gene.
	mark/deg_grem : path to MRF parameter files, usually a JLD/CSV file.
	mark/deg_hmm : path to HMM files, usually .hmm files.
	pop_size : size of population, i.e. number of individual HMM solutions to greedily optimize.
	frame : p1/p2.
	rel_change_thr : minimum threshold for the relative number of variants changing, used for setting a dynamic limit on the number of iterations.
	unchanged_thr : threshold for number of contiguous iterations without change to drop a seq from further optimization 
	max_iter : hard limit to the number of iterations during the MRF-based optimization
	mut_num : number of mutations per iteration as specified by CAMEOX input argument
	mut_len : length of mutations per iteration as specified by CAMEOX input argument	
    host_tid : NCBI Taxonomic ID for the host of the entanglement, used by the host generalization subsystem (default of 562 for E. coli)
    pll_weights:  control the selection of PLL weights for the entanglement pair with the following choices:
                equal (CAMEOS default), rand (paper and CAMEOX default), close2mark, close2deg. 
	X_range/Y_range: positions along gene to consider for double-encoding... useful if you want to restrict where double-encoding can occur.
						Ranges are defined in terms of the end of the sequence, so if you have a 70 aa sequence you want to start at position 80,
						set the range to be something like 150:151. Generally some buffer around this value (i.e. 150:160) is useful.
						Don't set, or set to false if you don't care where the genes overlap.
	actually_mrf: flag controlling whether doing MRF-based optimization
	debug: debugging argument
"""
function set_up_and_optimize(
    log_io, rand_barcode, paths, mark_name, deg_name, mark_grem, deg_grem,
    mark_hmm, deg_hmm, pop_size, frame, rel_change_thr, unchanged_thr, max_iter,
	mut_num, mut_len;
	host_tid = 562, pll_weights = "rand", gc_iter = 1,
	X_range = false, Y_range = false, actually_mrf = true, debug = 0,
	)

	@debug("Beginning run.")
	@debug(Libc.strftime(time()))
	@debug(Libc.strftime(time()))

	println("IMPORTANT: The random barcode on this run is: $rand_barcode")
	@debug("The random barcode on this run is: $rand_barcode")
	@debug("The host taxid for this run is: $host_tid")

	do_cull = false #culling reduces number of sequences we optimize over time. Just a trick to save time if you want top-performers only.
	full_even_window = 20
	half_window = div(full_even_window, 2) #window statistics to see if scores still decreasing (sometimes used).
	last_few_sf = Float64[]
	last_few_cs = Float64[]

	#Looks up mean/std. dev of family pseudolikelihoods, pre-computed.
	mu_mark, sig_mark = lookup.mu_sig(paths, deg_name)
	mu_deg, sig_deg = lookup.mu_sig(paths, mark_name)
	@debug("APLL stats: mu_mark=$(@sprintf("%.2f", mu_mark)), " * 
	 "sig_mark=$(@sprintf("%.2f", sig_mark)); mu_deg=$(@sprintf("%.2f", mu_deg)), " *
	 "sig_deg=$(@sprintf("%.2f", sig_deg))")

	# TODO: Expose gen_samples and candidate JLD file as input arguments
	gen_samples = true #this is true if starting from scratch (general case), false if we have some candidates to optimize ahead of time (untested CAMEOS legacy).
	population = types.Chromosome[]
	local mark_grem_prot, deg_grem_prot
	if !gen_samples
		population = load(joinpath(paths.input, "trpE_population.jld"))["population"][1:100]
		mark_gremodel, deg_gremodel = std_setup.short_set_up(mark_grem, deg_grem)
		deg_nNodes, mark_nNodes = deg_gremodel.nNodes, mark_gremodel.nNodes
	else
		#Generate population of sampled hmm starting points.
		if X_range == false || Y_range == false #we require both to be there...
			@debug("Doing standard full set up...")
			(mark_gremodel, deg_gremodel, population, mark_grem_prot,
			deg_grem_prot, real_frame) = std_setup.full_set_up(
				paths, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm,
				deg_grem, pop_size, rand_barcode, frame, host_tid)
		else
			@debug("The x range is $X_range")
			@debug("The y range is $Y_range")
			(mark_gremodel, deg_gremodel, population, mark_grem_prot,
			deg_grem_prot, real_frame) = std_setup.full_set_up(
				paths.output, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm,
				deg_grem, pop_size, rand_barcode, frame, host_tid;
				X_range, Y_range)
		end
		deg_nNodes, mark_nNodes = deg_gremodel.nNodes, mark_gremodel.nNodes
	end

	# Get (anti)pseudolikelihood for WT sequences
	function get_apll_wt()
		prot_seqs = bio_seq.load_fasta(joinpath(paths.input, "proteins.fasta"))
		local deg_seq, mark_seq
		# For mark
		mark_seq = prot_seqs[mark_name]
		@debug("$mark_name seq: $mark_seq")
		mark_ali_seq = uppercase(std_setup.align_consensus(mark_seq, mark_hmm))
		@debug("$mark_name aligned seq: $mark_ali_seq")            
		mark_wt_apll = mrf.psl(strip(mark_ali_seq, '*'),
		 mark_gremodel.w1, mark_gremodel.w2, true)
		@debug("$mark_name APLL: $mark_wt_apll")                       
		# For deg
		deg_seq = prot_seqs[deg_name]
		@debug("$deg_name seq: $deg_seq")
		deg_ali_seq = uppercase(std_setup.align_consensus(deg_seq, deg_hmm))
		@debug("$deg_name aligned seq: $deg_ali_seq")              
		deg_wt_apll = mrf.psl(strip(deg_ali_seq, '*'),
		 deg_gremodel.w1, deg_gremodel.w2, true)
		@debug("$deg_name APLL: $deg_wt_apll")
		# return WT APLLs
		return mark_wt_apll, deg_wt_apll     
	end
	mark_wt_apll::Float32, deg_wt_apll::Float32 = get_apll_wt()

	#Look through generated candidates, check their fitness values and their correctness.
	success_chrom, fdp, fmp, fitness_values, indie_deg_maps, indie_mark_maps, deg_ull, deg_pv_w1, deg_pv_w2, deg_prot_mat, mark_ull, mark_pv_w1, mark_pv_w2, mark_prot_mat = optimize.assess_founders(population, mark_gremodel.nNodes, mark_gremodel.w1, mark_gremodel.w2, deg_gremodel.nNodes, deg_gremodel.w1, deg_gremodel.w2)

	@debug("Num of successful seeds... $(length(success_chrom))")

	if actually_mrf
		founding_fitness_values = fitness_values[1:end]

		mark_base_energy = mrf.basic_energy_calc(
			mark_grem_prot, mark_gremodel.w1, mark_gremodel.w2)
		deg_base_energy = mrf.basic_energy_calc(
			deg_grem_prot, deg_gremodel.w1, deg_gremodel.w2)

		# Get the right PLL static weighting (rand is inside the loop later)
		static_1st_weight = 0.5  # case pll_weights = "equal"
		if pll_weights == "close2mark"
			static_1st_weight = 1.0
		elseif pll_weights == "close2deg"
			static_1st_weight = 0.0
		end

		# Now we load this info into an ExChrome type object
		cur_pop = types.ExChrome[]
		for i in eachindex(success_chrom)
			old = success_chrom[i]
			new_chrom = types.ExChrome(
				old.path, old.full_sequence, old.deg_nuc, indie_deg_maps[i],
				old.deg_trns, old.deg_trne, old.deg_d, old.deg_skip,
				old.deg_insert, old.mark_nuc, indie_mark_maps[i],
				old.mark_trns, old.mark_trne, old.mark_d, old.mark_skip,
				old.mark_insert, fdp[i], deg_base_energy,
				deg_prot_mat[i:i, 1:end],
				utils.aa_vec_to_seq(deg_prot_mat[i:i, 1:end]),
				deg_ull[i], deg_pv_w1[i], deg_pv_w2[i:i, 1:end], fmp[i],
				mark_base_energy, mark_prot_mat[i:i, 1:end],
				utils.aa_vec_to_seq(mark_prot_mat[i:i, 1:end]), mark_ull[i],
				mark_pv_w1[i], mark_pv_w2[i:i, 1:end],
				pll_weights == "rand" ? rand() : static_1st_weight,
				0) #So now everything is in an extended chromosome object.
		push!(cur_pop, new_chrom)
		end

		if debug > 1  # Save population from HMM but before MRF optimization
			outp = joinpath(paths.entangle, "saved_hmm_$(rand_barcode).jld")
			@debug("Let's save our HMM non-MRF-optimized population to ", outp, "...")
			saved_pre = population_to_saveable(
				cur_pop;
				mark_wt_apll=mark_wt_apll, deg_wt_apll=deg_wt_apll)
			FileIO.save(outp, "variants", saved_pre)
			println("HMM non-MRF-optimized variants saved to ", outp)
		end

		#There are occasionally differences between HMM/MRF models of the proteins, mostly because HMMs are more flexible with insertions.
		#This generates mappings between positions in either data type.
		mark_hmm_to_grem_map = std_setup.get_explicit_mapping(mark_grem_prot, mark_hmm)
		deg_hmm_to_grem_map = std_setup.get_explicit_mapping(deg_grem_prot, deg_hmm)

		@debug("Hey the explicit mapping for mark_hmm_to_grem is $mark_hmm_to_grem_map")
		@debug("Hey the explicit mapping for deg_hmm_to_grem is $deg_hmm_to_grem_map")

		# Set up some variables needed in the dynamic condition for the main CAMEOX loop
		changed_seq = length(cur_pop)
		rel_changed_seq = 1.0
		iter = 0;

		#Let's keep track of optimization history.
		history_matrix = zeros(Float32, round(Int64, pop_size * 1.2), max_iter)

		println("Beginning MRF-based optimization for $(length(cur_pop)) seeds...")
		flush(stdout)
		@debug("About to start optimization.")
		@debug(Libc.strftime(time()))

		stop_early = false #a condition to evaluate every (few) iteration(s) to see if it's worth continuing a run. When true, we clean up and move to next gene pair.
		not_worth_continuing = false

		#Okay let's set up the energy normal distributions...

		deg_energy_mu, deg_energy_sig = lookup.get_energy_params(paths, deg_name)
		mark_energy_mu, mark_energy_sig = lookup.get_energy_params(paths, mark_name)

		mark_energy_normal = Normal(mark_energy_mu, mark_energy_sig)
		deg_energy_normal = Normal(deg_energy_mu, deg_energy_sig)
		
		done_pop = types.ExChrome[]

		### MAIN LOOP: MRF-based optimization ###	
		if gc_iter > 0
			GC.gc(true) # Force a full GC before the loop
			GC.enable(false)
		end
		while (rel_changed_seq > rel_change_thr && iter < max_iter) #&& !(stop_early) # && minimum((fitness_values)) > 400)
			sfv = sort(fitness_values)
			#best_50 = sfv[50]

			if debug > 0 || iter % 25 == 0
				println("INFO: Step $(iter): $(rel_changed_seq*100)% of changed seqs [threshold = $(rel_change_thr*100)%]")
				flush(stdout)
				@debug("Iteration $iter: $(rel_changed_seq*100)% ($changed_seq) of changed seqs [threshold = $(rel_change_thr*100)%]")
				@debug("  Fitness stats: mean=$(mean((fitness_values))) std=$(std((fitness_values)))")
				@debug("  Best five were $(sfv[1:5])")
				@debug("  ... and worse five were $(sfv[end-5:end]).")
				flush(log_io)
			end

			#iterated conditional modes, multi-site keep-track = icm_multi_kt. keep track = store values of changes to sequence.
			if do_cull
				if iter < div(max_iter, 4)
					@debug("Cull for iter $iter : no cutoff")
					cur_pop, changed_seq = optimize.icm_multi_kt(
						cur_pop, deg_gremodel, mark_gremodel,
						mut_num = mut_num, mut_len = mut_len) #, deg_energy_normal, mark_energy_normal)
				elseif iter >= div(max_iter, 4) && iter < div(max_iter, 2)
					@debug("Cull for iter $iter : cutoff is $(sfv[div(length(sfv), 2)])")
					cur_pop, changed_seq = optimize.icm_multi_kt(
						cur_pop, deg_gremodel, mark_gremodel;
						score_cutoff = sfv[div(length(sfv), 2)],
						mut_num = mut_num, mut_len = mut_len)
				elseif iter >= div(max_iter, 2)
					@debug("Cull for iter $iter : cutoff is $(sfv[div(length(sfv), 4)])")
					cur_pop, changed_seq = optimize.icm_multi_kt(
						cur_pop, deg_gremodel, mark_gremodel;
						score_cutoff = sfv[div(length(sfv), 4)],
						mut_num = mut_num, mut_len = mut_len)
				end
			else
				cur_pop, changed_seq = optimize.icm_multi_kt(
					cur_pop, deg_gremodel, mark_gremodel;
					mut_num = mut_num, mut_len = mut_len)
			end
			if (gc_iter > 0) && (iter % gc_iter == 0)
				GC.enable(true)
				GC.gc(false) # Force GC after icm_multi_kt
				GC.enable(false)
			end

			fitness_values = optimize.assess_pop([done_pop; cur_pop])
			summed_fitness = sum(fitness_values) #one of the things we should be watching.

			push!(last_few_sf, summed_fitness)
			push!(last_few_cs, changed_seq)

			tmp_done_pop = types.ExChrome[]
			tmp_done_pop = filter(v->v.unchanged>unchanged_thr, cur_pop)
			append!(done_pop, tmp_done_pop)
			filter!(v->v.unchanged<=unchanged_thr, cur_pop)
			
			indiv_count = 1
			for indiv in fitness_values #history matrix stores values of function being optimized over time.
				history_matrix[indiv_count, iter + 1] = indiv
				indiv_count += 1
			end

			#This was originally just code to test the pseudolikelihoods decreasing correctly. It does.
			#Now it's still sorta somewhat useful to log a few values of sequences in both deg/mark as optimization goes.
			if (mod(iter, 50) == 0 && iter != 0)
				@debug("Second opinion, these should be ~0.")
				for i_count in 1:min(5, length(cur_pop))
					#for the record this is not a good way to generate random numbers.
					#doesn't go from 1 to end.
					i = round(Int64, rand() * (length(cur_pop) - 1) + 1) # it was just = i_count

					bel_deg = cur_pop[i].deg_prob
					bel_mrk = cur_pop[i].mark_prob
					#Sort of assuming that end-1 is because of a "*" at stop codon being added.
					if (length(strip(cur_pop[i].deg_seq, '*')) == deg_nNodes 
							&& length(strip(cur_pop[i].mark_seq, '*')) == mark_nNodes)
						cal_deg = mrf.psl(
							strip(cur_pop[i].deg_seq, '*'), deg_gremodel.w1,
							 deg_gremodel.w2, false)
						cal_mrk = mrf.psl(
							strip(cur_pop[i].mark_seq, '*'), mark_gremodel.w1,
							 mark_gremodel.w2, false)
						@debug("DEG: $bel_deg vs. $cal_deg")
						@debug("MRK: $bel_mrk vs. $cal_mrk")
						@debug(abs(bel_deg - cal_deg) < 1e-2)
						@debug(abs(bel_mrk - cal_mrk) < 1e-2)
					else
						@debug("Proteins not correct length!")
						@debug(cur_pop[i].mark_seq)
						@debug(cur_pop[i].deg_seq)
					end
				end
			end
			rel_changed_seq = changed_seq/length(success_chrom)
			iter += 1
			@debug("\n")
		end
		(gc_iter > 0) && GC.enable(true)

		#So we're done the iterated conditional modes step. Now we report everything...
		println("INFO: Optimization completed after iter $iter with $(rel_changed_seq*100)% changed seqs")
		flush(stdout)
		@debug("Done while loop after iter $iter with $(rel_changed_seq*100)% changed seqs!")
		@debug(Libc.strftime(time()))

		# Retrieve the remaining variants (if any) what were set to keep processing
		append!(done_pop, cur_pop)
		cur_pop = done_pop
		
		# Recalculate fitness for all the population and show stats
		fitness_values = optimize.assess_pop(cur_pop)
		sfv = sort(fitness_values)
		@debug("FINAL iteration $iter: $(rel_changed_seq*100)% ($changed_seq) of changed seqs [threshold = $(rel_change_thr*100)%]")
		@debug("  Fitness stats: mean=$(mean((fitness_values))) std=$(std((fitness_values)))")
		@debug("  Best five are $(sfv[1:5])")
		@debug("  ... and worse five are $(sfv[end-5:end]).")

		# History matrix
		@debug("Let us save the history matrix!")
		FileIO.save(
			joinpath(paths.entangle, "opt_his_mat_$(rand_barcode).jld"),
			"hist_mat", history_matrix)

		deg_significance = false
		mark_significance = false #this sees if we get any hits worth keeping at all. If yes we save full population.
		all_cal_deg_scores = Float64[]
		all_cal_mark_scores = Float64[]
		#If no we save just a few (top 10) individuals. These are top-6 overall + top-2 in deg and top-2 in mark.
		out_file = open(
			joinpath(paths.entangle, "all_final_fitness_$(rand_barcode).txt"), "w")
		write(out_file, "Ind. #\tMark Score\tDeg Score\tMark Sign\tDeg Sign\n")
		for last_indi in eachindex(cur_pop)
			rep_deg = cur_pop[last_indi].deg_prob
			rep_mrk = cur_pop[last_indi].mark_prob
			if (length(strip(cur_pop[last_indi].deg_seq, '*')) == deg_nNodes
					&& length(strip(cur_pop[last_indi].mark_seq, '*')) == mark_nNodes)
				cal_deg = abs(
					mrf.psl(strip(cur_pop[last_indi].deg_seq, '*'),
					deg_gremodel.w1, deg_gremodel.w2, false))
				cal_mrk = abs(
					mrf.psl(strip(cur_pop[last_indi].mark_seq, '*'),
					mark_gremodel.w1, mark_gremodel.w2, false))
				special_deg = "-"
				if cal_deg < mu_deg + sig_deg
					special_deg = "***!"
					deg_significance = true
				elseif cal_deg < mu_deg + 2*sig_deg
					special_deg = "**."
					deg_significance = true
				elseif cal_deg < mu_deg + 3*sig_deg
					special_deg = "*?"
					deg_significance = true
				elseif cal_deg < mu_deg + 4*sig_deg #we don't mind keeping these.
					deg_significance = true
				end
				special_mark = "-"
				if cal_mrk < mu_mark + sig_mark
					special_mark = "***!"
					mark_significance = true
				elseif cal_mrk < mu_mark + 2*sig_mark
					special_mark = "**."
					mark_significance = true
				elseif cal_mrk < mu_mark + 3*sig_mark
					special_mark = "*?"
					mark_significance = true
				elseif cal_mrk < mu_mark + 4*sig_mark #we don't mind keeping these.
					mark_significance = true
				end
				push!(all_cal_mark_scores, cal_mrk)
				push!(all_cal_deg_scores, cal_deg)
				write(out_file, "$(last_indi)\t$(cal_mrk)\t$(cal_deg)\t$special_mark\t$special_deg\n")
			else
				push!(all_cal_mark_scores, Inf)
				push!(all_cal_deg_scores, Inf)
				write(out_file, "$(last_indi)\t-1.0\t-1.0\t/\t/\n")
			end
		end
		close(out_file)

		cal_deg_p = sortperm(all_cal_deg_scores)
		cal_mark_p = sortperm(all_cal_mark_scores)
		cal_both_p = sortperm(fitness_values)

		@debug(Libc.strftime(time()))
		if FINALLY_SAVE_ALL_POPULATION || debug > 0  #deg_significance && mark_significance
			@debug("Let's save our optimized population...")
			saved_pop = population_to_saveable(
				cur_pop;
				mark_wt_apll=mark_wt_apll, deg_wt_apll=deg_wt_apll)
			FileIO.save(
				joinpath(paths.entangle, "saved_pop_$(rand_barcode).jld"), 
				"variants", saved_pop)
		else #We want just a subset.
			@debug("We'll save 12 interesting members of the population...")
			selected_pop = ExChrome[]
			for ijk in 1:3
				push!(selected_pop, cur_pop[cal_deg_p[ijk]])
				push!(selected_pop, cur_pop[cal_mark_p[ijk]])
			end
			for ijk in 1:6
				push!(selected_pop, cur_pop[cal_both_p[ijk]])
			end
			FileIO.save(
				joinpath(paths.entangle, "top_pop_$(rand_barcode).jld"),
				"variants", selected_pop)
		end

		@debug("And we'll also output our top twelve sequences...")
		out_file = open(joinpath(paths.entangle, "top_twelve_$(rand_barcode).fa"), "w")

		for ijk in 1:3
			deg_ind = cal_deg_p[ijk]
			write(out_file, ">TOP_DEG_d (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
			write(out_file, "$(cur_pop[deg_ind].deg_seq)\n")
			write(out_file, ">TOP_DEG_m (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
			write(out_file, "$(cur_pop[deg_ind].mark_seq)\n")
		end
		for ijk in 1:3
			deg_ind = cal_mark_p[ijk]
			write(out_file, ">TOP_MARK_d (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
			write(out_file, "$(cur_pop[deg_ind].deg_seq)\n")
			write(out_file, ">TOP_MARK_m (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
			write(out_file, "$(cur_pop[deg_ind].mark_seq)\n")
		end
		for ijk in 1:6
			deg_ind = cal_both_p[ijk]
			write(out_file, ">TOP_GEN_d (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
			write(out_file, "$(cur_pop[deg_ind].deg_seq)\n")
			write(out_file, ">TOP_GEN_m (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
			write(out_file, "$(cur_pop[deg_ind].mark_seq)\n")
		end
		close(out_file)

		@debug("Finally, we'll save some metadata of the run...")
		metadata_filename = joinpath(paths.entangle, "CAMEOX_metadata.csv") 
		if !isfile(metadata_filename)
			out_file = open(metadata_filename, "w")  
			write(out_file,
			 "rand_barcode,mark_name,deg_name,pop_size,frame,rel_change_thr,pll_weights,mark_wt_apll,deg_wt_apll,rel_changed_seq,iters,max_iters,sfv_top,sfv_end,mean_fitness,std_fitness,unchanged_thr,mut_num,mut_len,datetime,elapsed,threads,gc_iter\n")
		else
			out_file = open(metadata_filename, "a")                
		end
		println(out_file,
		 rand_barcode, ',', mark_name, ',', deg_name, ',', pop_size, ',', real_frame, ',',
		 rel_change_thr, ',', pll_weights, ',', mark_wt_apll, ',', deg_wt_apll, ',',
		 rel_changed_seq, ',', iter, ',', max_iter, ',', sfv[1], ',', sfv[end], ',', 
		 mean(fitness_values), ',', std(fitness_values), ',', unchanged_thr, ',',
		 mut_num, ',', mut_len, ',', round(now(), Dates.Second), ',', 
		 Dates.format(Dates.Time(Dates.Nanosecond(now() - INI_TIME)), "HH:MM:SS"), ',',
		 Threads.nthreads(), ',', gc_iter
		 )
		close(out_file)
	end
	@debug("CAMEOX DONE!")
	@debug(Libc.strftime(time()))
end

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table! s begin
		"tasks"
			help = "positional arg, path to TSV file with tasks (entanglements) to process"
			default = "aroB_infA_pf5_uref100_debug.txt"
		"--maxiter", "-m"
			help = "hard limit to the number of iterations in the MRF-based optimization"
			arg_type = Int
			default = 5000
		"--unchanged", "-u"
			help = "contiguous iters without change to drop a seq from further optimization"
			arg_type = Int
			default = 1000			
		"--nomrf"
			help = "skip the MRF-based optimization"
			action = :store_true		
		"--mutnum", "-n"
			help = "integer number of mutations or 'rand' (random for each iteration)"
			range_tester = (x-> x=="rand" || 0 < parse(Int64,x))
			default = "1"
 		"--mutlen", "-l"
			help = "integer length of the mutation or 'rand' (random for each iteration)"
			range_tester = (x-> x=="rand" || 0 < parse(Int64,x)) 
			default = "3"
		"--gciter"
			help = "number of MRF-based optimization iterations per GC event (0 disables)"
			arg_type = Int
			default = 1
		"--basedir", "-b"
			help = "base directory for data (defaults to CAMEOS behaviour of code dir)"
			arg_type = AbstractString
			default = "."		
		"--debug", "-g"
			help = "activate debugging mode (repeat for further debugging output)"
			action = :count_invocations
		"--num"
			help = "[CAMEOS legacy] experimentally used for running multiple jobs at once"
			default = "0"
		"--threads"
			help = "[CAMEOS legacy] used for parallelizing each job within a node"
			default = "1"

	end

	return parse_args(s)
end

function run_file()
    println("=-= CAMEOX = CAMEOs eXtended =-= v1.0.1 - Aug 2023 =-= LLNL =-=")
	flush(stdout)

	parsed_args = parse_commandline()
	task_file = parsed_args["tasks"]
	max_iter = parsed_args["maxiter"]
	unchanged_thr = parsed_args["unchanged"]
	no_mrf = parsed_args["nomrf"]
	mut_num = parsed_args["mutnum"]
	mut_len = parsed_args["mutlen"]
	gc_iter = parsed_args["gciter"]
	debug = parsed_args["debug"]
	num = parsed_args["num"]
	threads = parsed_args["threads"]
	basedir = parsed_args["basedir"]

	if debug > 0
		println("These CAMEOX command arguments are in effect:")
		for key in keys(parsed_args)
			println("  > '", key, "' is ", parsed_args[key])
		end
		if debug > 2
			GC.enable_logging(true)
		end
	end

	gc_iter > 1 && println(
		"CAUTION: GC_iter=", gc_iter, " is larger than default. Check memory usage!")
	flush(stdout)

	RUN_I = parse(Int64, num)
	NUM_THREADS = parse(Int64, threads)

	if RUN_I > 1 && NUM_THREADS > 1
		sleep(rand() * 2) #randomly delay start to de-synchronize starts of runs.
	end

	if isfile(task_file)
		in_file = open(task_file)
		in_read = readlines(in_file)
		close(in_file)

		#This is an additional file that records errors independent of individual log files.
		problem_file = open(joinpath(basedir, "problem_runs_$(RUN_I).txt"), "a")

		line_count = 0
		for line in in_read #
			line_count += 1
			if line[1] != '#' && (line_count % NUM_THREADS == RUN_I)
                rand_barcode = Random.randstring()
				run_args = split(line, '\t')
                try
                    out_dir, short, long, short_jld, long_jld, short_hmm,
                    long_hmm, pop_size, frame, rel_change_thr, host_tid,
                    pll_weights = run_args

                    host_tid = parse(Int64, host_tid)
                    if host_tid == 0  #If default taxid, then use E. coli taxid
                        host_tid = 562
                    end

                    if !(pll_weights in [
                        "equal", "rand", "close2mark", "close2deg"])
                        throw(DomainError(
                            pll_weights,
                            "this PLL optimization choice is unknown!"))
                    end
                    
					out_path = joinpath(basedir, out_dir)
					entangle_path = joinpath(out_path, "$(short)_$(long)_$frame")
					if !(isdir(entangle_path))
						run(`mkdir -p $entangle_path`)
					end
					debug > 0 && println("DEBUG: Entangle output path is $entangle_path")
					# Define important paths in a Path instance
					paths = types.Paths((
						base=basedir,
						input=basedir,
						energies=joinpath(basedir, "energies/"),
						psls=joinpath(basedir, "psls/"),
						output=out_path,
						entangle=entangle_path
						))				

					log_path = joinpath(entangle_path, "log_$(rand_barcode).txt")
					log_io = open(log_path, "w+")
					debug > 0 && println("DEBUG: Full log stored on $log_path")

					logger = SimpleLogger(log_io, Logging.Debug)
					global_logger(logger)

					with_logger(logger) do
						pop_size = parse(Int64, pop_size)
						rel_change_thr = parse(Float64, rel_change_thr)
						set_up_and_optimize(log_io, rand_barcode, paths,
                                            short, long,
                                            joinpath(paths.input, short_jld),
											joinpath(paths.input, long_jld),
											joinpath(paths.input, short_hmm),
											joinpath(paths.input, long_hmm),
                                            pop_size, frame,
                                            rel_change_thr, unchanged_thr, max_iter,
											mut_num, mut_len;
											host_tid = host_tid, pll_weights = pll_weights,
											gc_iter = gc_iter,
											actually_mrf = !no_mrf, debug = debug,
											)
					end

					flush(log_io)
					close(log_io)
				catch y
					write(problem_file,
                          "BC $rand_barcode ==> Problem $y processing line: $line\n")
                    flush(problem_file)
					if debug > 0
						rethrow()
					end
				end
			end
		end
		close(problem_file)
	else
		println("$task_file does not exist. Exiting.")
	end

end

const INFO_ITER_STEP = 25 # Iterations per info cycle (not debugging)
const FINALLY_SAVE_ALL_POPULATION = true # Save all the population or just ~10 best variants
const INI_TIME = now() # Initial time for elapsed time calculations
@time run_file()

#
# Profiling (uncomment lines below if you are interested in profiling the code)
#
#@profile run_file()
#open("CAMEOX_prof.txt", "w") do s
#    Profile.print(IOContext(s, :displaysize => (24, 500)))
#end
#win = Gtk.Window("gtkwait")
#ProfileView.view()
#if !isinteractive()
#    @async Gtk.gtk_main()
#    Gtk.waitforsignal(win,:destroy)
#end

end
