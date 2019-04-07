import "deepvariant.wdl" as dvtool

workflow DV {
	
	File input_bam
#	File output_prefix
	File Ref
	File bedFile

	call dvtool.make_examples {
		input: input_bam = input_bam, ref_fasta = Ref, BedFile = bedFile, Examples = "examples",
	}
	call dvtool.call_variants {
		input: Examples = make_examples.ExamplesOutput, name_flag= "call_variants"
	}
	call dvtool.post_process {
		input: ref_fasta = Ref, InputFile = call_variants.out, FinalOutput = "deepVariant_raw",
	}
	output {
		File vcf = post_process.Output
	}
}
