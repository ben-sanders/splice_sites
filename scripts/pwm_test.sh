# 6. PWM
# ------

# Check against the derived PWMs from the S&S analysis subsection
# input format is very similary to MESpy.

cd ../refseq_sites

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="REFSEQ";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 50, 15)"\t"substr($4, 50, 15)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 38, 15)"\t"substr($4, 38, 15)"\t"$5"\t"$6 }}' \
> input.refseq_sites.pwm

python ../scripts/python/pwm_score.py input.refseq_sites.pwm > out.refseq_sites.pwm

# tidy up
mv input.refseq_sites.pwm complete/

cd ../random_sites

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="RANDOM";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 50, 15)"\t"substr($4, 50, 15)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 38, 15)"\t"substr($4, 38, 15)"\t"$5"\t"$6 }}' \
> input.random_sites.pwm

python ../scripts/python/pwm_score.py input.random_sites.pwm > out.random_sites.pwm

# tidy up
mv input.random_sites.pwm complete/
