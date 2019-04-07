task strelkaGermline {
	File normal_bam
	File Ref_file
	String outprefix
	String tmpDIR = "strelkaTMP_" + outprefix 
	File calling_intervals

	command <<<
		/opt/bin/configureStrelkaGermlineWorkflow.py \
		--bam ${normal_bam} \
		--ref ${Ref_file} \
		--callRegions ${calling_intervals}
		--runDir ${tmpDIR}

		${tmpDIR}/runWorkflow.py -m local -j 10
		mv ${tmpDIR}/results/variants/variants.vcf.gz ${outprefix}.strelka.germline.vcf.gz
		mv ${tmpDIR}/results/variants/genome.S1.vcf.gz ${outprefix}.strelka.genome.germline.vcf.gz
	>>>
	runtime {
		docker: "strelka:v1"
		cpu:"2"
		memorty:"20G"
	}

	output {
		File strelkaGermlineVCF = "${outprefix}.strelka.germline.vcf.gz"
	}
}

