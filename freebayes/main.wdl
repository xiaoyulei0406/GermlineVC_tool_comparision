#import "pre_processing.wdl" as maptasks
import "freebayes.wdl" as freebayestasks
workflow freebayes {
	
	File input_bam
	String output_prefix
	File Ref
	File calling_intervals
	call freebayestasks.Germlinecall as Germlinecall {
		input: bam=input_bam, Ref_file = Ref, outprefix = output_prefix, bedFile = calling_intervals
	}

}
