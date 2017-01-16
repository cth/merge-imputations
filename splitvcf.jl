#!/home/fng514/bin/julia
#$ -S /home/fng514/bin/julia
#$ -cwd 
# Christian Theil Have, 2017.
#
using ArgParse
using BGZFStreams

@everywhere function process_line(line1,lineno) 
	if (ismatch(r"^#", line1)) 
		(0, line1)
	else
		(lineno, line1)
	end
end

# FIXME: We may need to do this with gzip readers instead
function main(args)
	s = ArgParseSettings("Example 2 for merge.jl: " *  # description
		"flags, options help, " *
		"required arguments.")

	@add_arg_table s begin
		"input-vcf"
		"--lines"
		arg_type = Int
		default = 1000
		"--prefix"
	end

	parsed_args = parse_args(s) # the result is a Dict{String,Any}
	println("Parsed args:")
	for (key,val) in parsed_args
		println("  $key  =>  $(repr(val))")
	end

	vcf1=BGZFStream(parsed_args["input-vcf"])
	out = BGZFStream(parsed_args["prefix"],"w")
	lines1 = eachline(vcf1)

	state1 = start(lines1)

	keep_positions1 = nothing

	header_lines = []
	content_lines = []

	fileno = 1
	lineno = 0

	while !done(lines1, nothing)
		(current_line1, state1) = next(lines1, state1)
		(lineno, processed_line)  =  process_line(current_line1, lineno+1) 

		if lineno == 0
			push!(header_lines, processed_line)
		elseif lineno % parsed_args["lines"] == 0
			push!(content_lines,processed_line)
			out_name = string(parsed_args["prefix"],lpad(fileno, 5, 0),".vcf.gz")
			out = BGZFStream(out_name, "w")
			for line in header_lines
				write(out,line)
			end
			for line in content_lines
				write(out,line)
			end
			close(out)
			content_lines = []
			fileno = fileno + 1
			println("wrote $out_name")
		else
			push!(content_lines,processed_line)
		end
	end

	if length(content_lines) > 0
		out_name = string(parsed_args["prefix"],lpad(fileno, 5, 0),".vcf.gz")
		out = BGZFStream(out_name, "w")
		for line in header_lines
			write(out,line)
		end
		for line in content_lines
			write(out,line)
		end
		close(out)
		println("wrote $out_name")
	end

	close(vcf1)
end

main(ARGS)
