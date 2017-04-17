dir=$1

echo "dir ${dir}"

for member in ${dir}/*; do
    if [[ ${member} == "${dir}/Sample"* ]]; then
        echo "accessing ${member}"
        analysis_dir="${member}/analysis"
        for vcf_file in ${analysis_dir}/*; do
            if [[ ${vcf_file} == *"filtered.vcf.gz" ]]; then
                echo "removing ${vcf_file}"
                rm -rf ${vcf_file}
            fi
            if [[ ${vcf_file} == *".vcf.gz.tbi" ]]; then
                echo "removing ${vcf_file}"
                rm -rf ${vcf_file}
            fi
            if [[ ${vcf_file} == *"Calls.vcf.gz" ]]; then
                echo "processing ${vcf_file}"
                vcf_filename="${vcf_file:0:(${#vcf_file}-3)}"
                echo "vcf_filename = ${vcf_filename}"
                echo "zcating ${vcf_file} into ${vcf_filename}"
                zcat ${vcf_file} > ${vcf_filename} 
                echo "removing ${vcf_file}"
                rm -rf ${vcf_file}
                echo "bgzipp ${vcf_filename}"
                bgzip ${vcf_filename}
                echo "tabix ${vcf_file}"
                args="-p vcf"
                tabix ${args} ${vcf_file}
                echo "tabix file created ${vcf_file}"
            fi
        done
    fi
done
