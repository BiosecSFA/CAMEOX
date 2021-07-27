"""
host: struct dealing with different hosts (E.coli, pf5, and more)
"""

module host

using DataFrames, CSV

aa2codons = Dict{Char, Array{String, 1}}(
    'A' => ["GCT", "GCC", "GCA", "GCG"],
    'R' => ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    'N' => ["AAT", "AAC"],
    'D' => ["GAT", "GAC"],
    'B' => ["AAT", "AAC", "GAT", "GAC"],  # N or D
    'C' => ["TGT", "TGC"],
    'Q' => ["CAA", "CAG"],
    'E' => ["GAA", "GAG"],
    'Z' => ["CAA", "CAG", "GAA", "GAG"],  # Q or E
    'G' => ["GGT", "GGC", "GGA", "GGG"],
    'H' => ["CAT", "CAC"],
    'I' => ["ATT", "ATC", "ATA"],
    'L' => ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
    'K' => ["AAA", "AAG"],
    'M' => ["ATG"],
    'F' => ["TTT", "TTC"],
    'P' => ["CCT", "CCC", "CCA", "CCG"],
    'S' => ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    'T' => ["ACT", "ACC", "ACA", "ACG"],
    'W' => ["TGG"],
    'Y' => ["TAT", "TAC"],
    'V' => ["GTT", "GTC", "GTA", "GTG"],
    '*' => ["TAA", "TGA", "TAG"]
)

cut = DataFrame()  # Global to be able to use eval in struct Host

struct Host{T<:Integer}
    taxid::T
    name::String
    opt_codon::Dict{Char, String}

    function Host{T}(tid::T) where {T<:Integer}
        # Select supported host
        name::String = ""
        if tid == 562
            name = "Escherichia coli"
        elseif tid == 220664 || tid == 1218948  # Heterotypic synonyms
            tid = 220664
            name = "Pseudomonas sp. Pf-5"
        else
            name = "Unknown taxon"
        end

        # Parse CUT from file to get the aa => codon optimization dictionary 
        aa2opt_codon = Dict{Char, String}()
        filenamecut::String = "CUT_$(tid).tsv"
        if !isfile(filenamecut)
            error("Missing Codon Usage Table file $filenamecut for taxid $(tid)!")
        end
        global cut = CSV.read("CUT_$(tid).tsv", DataFrame)  # CUT = Codon Usage Table
        for (aa, codons) in aa2codons
            # Build logical expression for each AA attending to codon redundancy
            #  e.g. for F/Phe: (dat.codon .== "TTT") .| (dat.codon .== "TTC")
            logic_str::String = ""
            for (cnt, codon) in enumerate(aa2codons[aa])
                if cnt > 1
                    logic_str = logic_str * " .| " 
                end
                logic_str = logic_str * "(cut.codon .== \"$codon\")"
            end
            # Use metaprogramming to parse and eval the logical expression
            cut_aa = cut[eval(Meta.parse(logic_str)),:]
            mxval, mxidx = findmax(cut_aa.freq) # Get index of rel freq max
            opt_codon = cut_aa[mxidx,1]  # Get optimum codon using index
            #println(aa, ' ', cut_aa, ' ', opt_codon)
            aa2opt_codon[aa] = opt_codon
        end       
        println("INFO: Host: Parsed Codon Usage Table for taxid $tid ($name)")

        new(tid, name, aa2opt_codon)  # Build new object of type Host
    end
end


function optimize_codons(host, aa_seq::AbstractString)
	return join([host.opt_codon[aa] for aa in aa_seq])
end

end


