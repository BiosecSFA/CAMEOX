"""
outparser: CAMEOX output parser
"""

using JLD # Should be placed here too to avoid later crash loading JLD file

module main

include("optimize.jl")
include("types.jl") 

using ArgParse, JLD

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table s begin
        "mark_gene"
		    help = "1st positional argument, name of gene 'mark'"
            arg_type = String
			default = "cysJ"
        "deg_gene"
		    help = "2nd positional argument, name of gene 'deg'"
            arg_type = String
			default = "infA"
		"runid"
			help = "3rd positional argument, run id code created by CAMEOX"
            arg_type = String
            default = "LgolUs7k"
        "--frame"
			help = "optional argument, CAMEOX frame"
            arg_type = String
		    default = "p1"     
        "--fasta"
            help = "generate also a FASTA file with specific content: all, fullseq, mark_gene, mark_seq, deg_gene, deg_seq"
            arg_type = String
	end
	return parse_args(s)
end

function outparse_cameos()
    println("=-= CAMEOX output parser =-= v0.6 - Jul 2021 =-= by LLNL =-=")

    # Parse arguments
   	parsed_args = parse_commandline()
	mark_gene = parsed_args["mark_gene"]
	deg_gene = parsed_args["deg_gene"]
	runid = parsed_args["runid"]
	frame = parsed_args["frame"]
    fasta = parsed_args["fasta"]

    # Get complete path for input and output files
    subdir = string(mark_gene, "_", deg_gene, "_", frame)
    jld_file = string("saved_pop_", runid, ".jld")
    csv_file = string("summary_", runid, ".csv")
    fa_file = string("variants_", runid, ".fasta")
    if fasta != nothing && fasta != "all"
        fa_file = string("variants_", runid, "_", fasta, ".fasta")
    end
    in_path = string("output/", subdir, "/", jld_file)
    out_path = string("output/", subdir, "/", csv_file)
    fa_path = string("output/", subdir, "/", fa_file)
    if fasta == "fullseq" # Rename with the long field
        fasta == "full_sequence"
    end

    # Load JLD file
    print("Loading data from JLD file ", in_path, " ... ")
    variants = load(in_path)["variants"];
    println("OK!")

    # Save CSV file
    print("Saving data to CSV file ", out_path, " : ")
    parsed = 0
    open(out_path, "w") do io
        println(io, string("full_seq,",
                           mark_gene,"_psls,",
                           mark_gene,"_seq,",
                           deg_gene,"_psls,",
                           deg_gene,"_seq"))
        for var in variants
            println(io, var.full_sequence,",",
                    var.mark_prob,",",
                    var.mark_seq,",", #var.mark_nuc,",",
                    var.deg_prob,",",
                    var.deg_seq) #,",",var.deg_nuc)
            parsed += 1
            if parsed % 200 == 0
                print(".")
                flush(stdout)
            end
        end
    end
    println(" OK!")
    println(parsed, " variants parsed for ", runid)

    # Save fasta file
    if fasta != nothing
        print("Saving ", fasta," to FASTA file ", fa_path, " : ")
        variant = 1
        open(fa_path, "w") do io

            function print2fasta(var, field::String, hdrverb::Bool=true,
                                 nucs::Bool=false)
                # Prints header and sequence to fasta file

                # Prints header
                if hdrverb
                    println(io, ">CAMEOX run:", runid,
                            " mark_gene:", mark_gene,
                            " deg_gene:", deg_gene,
                            " variant:", variant, "/", parsed,
                            " ", string(field))
                elseif occursin("mark", field)
                    println(io, ">", mark_gene, "_", variant, "_", runid)
                elseif occursin("deg", field)
                    println(io, ">", deg_gene, "_", variant, "_", runid)               
                else
                    println(io, ">", field, "_", variant, "_", runid)
                end
                
                # Prints sequence
                if nucs
                    println(io, getfield(getfield(var, Symbol(field)), :nucs))
                else
                    println(io, getfield(var, Symbol(field)))
                end
                return nothing
            end

            # Loop over variants printing depending on selection
            for var in variants
                if fasta == "full_sequence"
                    print2fasta(var, "full_sequence", true)  # equal to mark_nuc.nucs
                elseif fasta == "all" || fasta == "mark_seq"
                    print2fasta(var, "mark_seq", (fasta == "all"))
                end
                if fasta == "all" || fasta == "mark_nuc"
                    print2fasta(var, "mark_nuc", (fasta == "all"), true)
                end
                if fasta == "all" || fasta == "deg_seq"                    
                    print2fasta(var, "deg_seq", (fasta == "all"))
                end
                if fasta == "all" || fasta == "deg_nuc"                    
                    print2fasta(var, "deg_nuc", (fasta == "all"), true)
                end
                variant += 1
                if variant % 200 == 0
                    print(".")
                end
            end
        end
        println(" OK!")
        println(variant-1, " variants in FASTA file ", fa_file)    
    end
end
    
@time outparse_cameos()
end
