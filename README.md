# genome

usage: main.py [-h] [-d DIRECTORY] [-f] [-m] [-o OUTPUT] [-ht]

optional arguments: \n
  -h, --help            show this help message and exit \n
  -d DIRECTORY, --directory DIRECTORY, --dir DIRECTORY \n
                        parent directory of the samples i.e. Sample_054 \n
  -f, --filter          use argument to keep the following columns: CHROM,\n
                        POS, REF, ALT, genotype info\n
  -m, --merge           use argument to merge all filtered.vcf.gz files in the \n
                        parent directory \n
  -o OUTPUT, --output OUTPUT\n
                        directory where to output the filtered or merged files\n
  -ht, --homozygous_test\n
                        use argument to collect homozygote statistics\n

