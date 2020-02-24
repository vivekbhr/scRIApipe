#!/bin/bash
# input file:
#     - $1 = input = "transcripts_quant/{sample}/output.correct.sort.txt"
# output files:
#     - $2 = mtx = "transcripts_quant/{sample}/eq_counts/output.mtx"
#     - $3 = ec = "transcripts_quant/{sample}/eq_counts/output.ec.txt"
#     - $4 = bc = "transcripts_quant/{sample}/eq_counts/output.barcodes.txt"
# temp files:
#     - $5 =  bc_index = temp("transcripts_quant/{sample}/eq_counts/bc_index.txt"),
#     - $6 = ec_index = temp("transcripts_quant/{sample}/eq_counts/ec_index.txt"),
#     - $7 = header = temp("transcripts_quant/{sample}/eq_counts/mtx_header.txt"),
#     - $8 = tmp1 = temp("transcripts_quant/{sample}/eq_counts/tmp1.txt"),
#     - $9 = tmp2 = temp("transcripts_quant/{sample}/eq_counts/tmp2.txt")

# count unique RNA fragments based on unique combination of Barcode, UMI and EC
awk -F "\t" 'BEGIN{OFS=FS} {print($1, $3)}' $1 | sort | uniq -c | \
awk -F " " 'BEGIN{OFS="\t"} {print($2, $3, $1)}' > $8

# change index of tmp.mtx from barcode to number
cut -f 1 $8 | sort | uniq | awk 'BEGIN{OFS="\t"} {print($1, NR)}' > $5
awk 'FNR==NR{a[$1]=$2; next} {print(a[$1], $2, $3)}' $5 $8 > $9
cut -f 1 $5 > $4

# change index of ec
cut -d " " -f 2 $9 | sort -n | uniq | awk 'BEGIN{OFS="\t"} {print($1, NR)}' > $6
awk 'FNR==NR{a[$1]=$2; next} {print($1, a[$2], $3)}' $6 $9 > $8

# make header for mtx file output
tail -n 1 $5 | cut -f 2 > $7
cut -d " " -f 2 $8 | sort -n -r | head -n 1 >> $7
awk 'END{{print(NR)}}' $8 >> $7
echo "%%MatrixMarket matrix coordinate real general" > $2
awk 'BEGIN{ORS=" "}1' $7 >> $2
echo "" >> $2
cat $8 >> $2

cut -f 1 $6 > $3

echo 'done'
