task bqsr {
	
	File bam
	File Ref_file
	File knowsite1
	File knowsite2
	File knowsite3
	String sample_name_flag

	command <<<
		gatk BaseRecalibrator --input ${bam} --output ${sample_name_flag}.bqsrtable --reference ${Ref_file} --known-sites ${knowsite1} --known-sites ${knowsite2} --known-sites ${knowsite3}

		gatk ApplyBQSR -R ${Ref_file} -I ${bam} -bqsr ${sample_name_flag}.bqsrtable -O ${sample_name_flag}.bqsr.bam
	>>>

	output {
		File recal_table = "${sample_name_flag}.bqsrtable"
		File recalibrated_bam = "${sample_name_flag}.bqsr.bam"
		File outbai = "${sample_name_flag}.bqsr.bam.bai"
	}
	
	runtime{
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"50G"
	}
}

task haplotypecaller {
	File Ref_file
	String sample_name_flag
	File recal_table
	File bam
	File dbSNP
#	String outdir
	File intervals

	command {
		gatk HaplotypeCaller -R ${Ref_file} \
		-I ${bam} \
		--dbsnp ${dbSNP} \
		-L ${intervals} \
		-O ${sample_name_flag + ".vcf"}
	}
	output { 
		File vcf = "${sample_name_flag}.vcf"
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"30G"
	}
}

task GenotypeGVCFs {
	File Ref_file
#	Array[File] GVCFs
	File my_database
#	String outdir
	String sample_name_flag

	command <<<
		tar -xf ${my_database}
		WORKSPACE=$( basename ${my_database} .tar)
		echo $WORKSPACE

		java -jar /gatk/gatk.jar GenotypeGVCFs -R ${Ref_file} \
		-V gendb://$WORKSPACE \
		-O ${sample_name_flag}.vcf
	>>>

	output { 
		File vcf = "${sample_name_flag}.vcf"
	}

	runtime {
		docker:"gatk4/wd:latest"
		cpu:"4"
		memory:"30G"
	}
}
