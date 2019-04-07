
task Germlinecall {
	File bam
	File Ref_file
	String outprefix
	String tmpDIR = "strelkaTMP_" + outprefix
	File bedFile

	command <<<
		/usr/local/bin/freebayes \
		-f ${Ref_file} \
		${bam} \
		--targets ${bedFile} > ${outprefix}.vcf
	>>>
	runtime {
		docker: "freebayes:v1"
		cpu:"2"
		memorty:"20G"
	}

	output {
		File GermlineVCF = "${outprefix}.vcf"
	}
}

