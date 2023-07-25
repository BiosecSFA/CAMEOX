"""
lookup : keeps track of code that's main function is to look up external data.
"""

module lookup

using DelimitedFiles, Statistics

function get_energy_params(paths, gene_name) #NOTE: extra dependency
	deg_energy_file = open(joinpath(paths.energies, "energy_$(gene_name).txt"))
	deg_energy_read = readlines(deg_energy_file)
	close(deg_energy_file)

	#energies = deg_energy[strip(split(der, '\t')[2]) for der in deg_energy_read]
	energies = map(xx -> parse(Float64, xx), deg_energy_read)
	mean_energy = mean(energies)
	std_energy = std(energies)
	return mean_energy, std_energy
end

function mu_sig(paths, gene_name) #NOTE: extra dependency
	raw_data = readdlm(joinpath(paths.psls, "psls_$(gene_name).txt"))
	raw_values = map(Float64, raw_data[1:end, 1])
	mu = mean(raw_values)
	sig = std(raw_values)
	return mu, sig
end

end
