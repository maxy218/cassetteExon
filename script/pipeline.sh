bin_dir=$(cd "$(dirname "$0")"; cd ../bin; pwd)
echo "bin_dir: " $bin_dir

conf_name=$1
gene_annotation=`grep "^gene_annotation" ${conf_name} | awk -F"=" '{print $2}'`
exon_annotation=`grep "^exon_annotation" ${conf_name} | awk -F"=" '{print $2}'`
expression_type=`grep "^expression_type" ${conf_name} | awk -F"=" '{print $2}'`
expression_file=`grep "^expression_file" ${conf_name} | awk -F"=" '{print $2}'`
reference_genome_dir=`grep "^reference_genome_dir" ${conf_name} | awk -F"=" '{print $2}'`
sf_motif_map=`grep "^sf_motif_map" ${conf_name} | awk -F"=" '{print $2}'`

log_file="log.txt"
error_info="ERROR! pipeline is stopped by some reason. "
error_info=${error_info}"Please refer to ${log_file} to find more information."

#if [ ! -f "${gene_annotation}" ]; then
#  echo "file ${gene_annotation} doesn't exist!"
#fi

#if [ ! -f "${exon_annotation}" ]; then
#  echo "file ${exon_annotation} doesn't exist!"
#fi

output_dir="./output_cassetteExon/"
if [ ! -d ${output_dir} ]; then
  mkdir ${output_dir}
fi

tmp_output_dir=${output_dir}"/tmp/"
if [ ! -d ${tmp_output_dir} ]; then
  mkdir ${tmp_output_dir}
fi

#execute the annotation preprocess
echo "preprocessing gene annotation..." | tee -a ${log_file}
echo ${bin_dir}"/anno_preproc" ${gene_annotation} ${tmp_output_dir}"/filtered_gtf_anno.gtf" ${tmp_output_dir}"/exon_info_file" >> ${log_file}
#`${bin_dir}"/anno_preproc" ${gene_annotation} ${tmp_output_dir}"/filtered_gtf_anno.gtf" ${tmp_output_dir}"/exon_info_file"` 2>&1 | tee -a ${log_file}
`${bin_dir}"/anno_preproc" ${gene_annotation} ${tmp_output_dir}"/filtered_gtf_anno.gtf" ${tmp_output_dir}"/exon_info_file"`
ret=$?
if [ ${ret} -ne 0 ]; then
  echo ${error_info}
  exit ${ret}
fi

#execute the 
echo "dealing with cassetteExon..." | tee -a ${log_file}
echo ${bin_dir}"/cassetteExon" ${tmp_output_dir}"/filtered_gtf_anno.gtf" ${exon_annotation} ${expression_type} ${expression_file} ${output_dir}"/exon_inclusion_level" ${tmp_output_dir}"/exon_info_file" ${tmp_output_dir}"/output_5_regions" >> ${log_file}
${bin_dir}"/cassetteExon" ${tmp_output_dir}"/filtered_gtf_anno.gtf" ${exon_annotation} ${expression_type} ${expression_file} ${output_dir}"/exon_inclusion_level" ${tmp_output_dir}"/exon_info_file" ${tmp_output_dir}"/output_5_regions"
ret=$?
if [ ${ret} -ne 0 ]; then
  echo ${error_info}
  exit ${ret}
fi

#get 5 regions' sequence
echo "geting the sequence of 5 regions(UU, UD, exon, DU, DD)..." | tee -a ${log_file}
echo ${bin_dir}"/get_seq" ${reference_genome_dir} ${tmp_output_dir}"/output_5_regions" ${tmp_output_dir}"/out_region_seq" >> ${log_file}
${bin_dir}"/get_seq" ${reference_genome_dir} ${tmp_output_dir}"/output_5_regions" ${tmp_output_dir}"/out_region_seq"  2>&1 | tee -a ${log_file}
ret=$?
if [ ${ret} -ne 0 ]; then
  echo ${error_info}
  exit ${ret}
fi

#get the motif hits
echo "getting the motif hits..." | tee -a ${log_file}
echo ${bin_dir}"/find_motif" ${sf_motif_map} ${tmp_output_dir}"/out_region_seq" ${tmp_output_dir}"/filtered_motif" ${output_dir}"/out_exon_motif_hit" ${output_dir}"/out_exon_sf_hit" >> ${log_file}
${bin_dir}"/find_motif" ${sf_motif_map} ${tmp_output_dir}"/out_region_seq" ${tmp_output_dir}"/filtered_motif" ${output_dir}"/out_exon_motif_hit" ${output_dir}"/out_exon_sf_hit" 2>&1 | tee -a ${log_file}
ret=$?
if [ ${ret} -ne 0 ]; then
  echo ${error_info}
  exit ${ret}
fi



