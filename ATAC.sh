### Task Trim Adapter
fastp -i *_1.fq.gz -I *_2.fq.gz -o *_1.trimmed.fq.gz -O *_2.trimmed.fq.gz --detect_adapter_for_pe --thread 16 --json *_fastp.json --html *_fastp.html 2 > *_fastp.log

fastqc -t 20 *_*.fq.gz -o */QC_result

###Task Align
bowtie2 -p 30 -X 2000 -x */Macaca_index -1 *_1.trimmed.fq.gz  -2 *_2.trimmed.fq.gz -S *_align.sam --met-file *_align_result.log &> *_align.log
