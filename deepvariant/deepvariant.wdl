
task make_examples {
	File BedFile
	File input_bam
	File ref_fasta
	String Examples

	command {
		python /home/bin/make_examples.zip \
		--mode calling \
		--ref ${ref_fasta} \
		--reads ${input_bam} \
		--examples ${Examples}.tfrecord.gz \
		--regions ${BedFile}
        }

	output {
		File ExamplesOutput = "${Examples}.tfrecord.gz"
	}
    runtime {
		docker: "dajunluo/deepvariant:latest"
		cpu:"20"
		memory:"30G"
    }
}

task call_variants {
	File Examples
	String name_flag

	command {
		python /home/bin/call_variants.zip \
		--outfile ${name_flag}.tfrecord.gz \
		--examples ${Examples} \
		--checkpoint /home/models/model.ckpt
	}

	output {
		File out = "${name_flag}.tfrecord.gz"
	}

	runtime {
		docker: "dajunluo/deepvariant:latest"
		cpu:"20"
		memory:"30G"
	}
}

task post_process {
	File ref_fasta
	File InputFile
	String FinalOutput

	command {
		python /home/bin/postprocess_variants.zip \
		--ref ${ref_fasta} \
		--infile ${InputFile} \
		--outfile ${FinalOutput}.vcf.gz
	}

	output {
		File Output = "${FinalOutput}.vcf.gz"
	}
    runtime {
		docker: "dajunluo/deepvariant:latest"
		cpu:"20"
		memory:"30G"
    }
}

task GatherVCFs {
	Array[File] Input_Vcfs
	String Output_Vcf_Name

	command {
		java -jar /gatk/gatk.jar \
		MergeVcfs \
		-I ${sep=" -I " Input_Vcfs} \
		-O ${Output_Vcf_Name}.vcf.gz \
		--CREATE_INDEX true
	}

	output {
		File output_vcfs = "${Output_Vcf_Name}.vcf.gz"
	}

    runtime {
		docker: "gatk4/wd:latest"
		cpu:"4"
		memory:"30G"
    }
}
