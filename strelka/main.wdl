import "strelka.wdl" as strelkatasks
workflow strelka {
	
	File input_bam
	String output_prefix
	File Ref
	File calling_intervals

	call strelkatasks.strelkaGermline as strelkaGermline {
		input: normal_bam=input_bam, Ref_file = Ref, outprefix = output_prefix, calling_intervals = calling_intervals
	}
	output {
		File vcf = strelkaGermline.strelkaGermlineVCF
	}
}
