# genome
```
usage: main.py [-h] [-d DIRECTORY] [-f] [-m] [-o OUTPUT] [-ht]

optional arguments: 
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY, --dir DIRECTORY 
                        parent directory of the samples i.e. Sample_054 
  -f, --filter          use argument to keep the following columns: CHROM,
                        POS, REF, ALT, genotype info
  -m, --merge           use argument to merge all filtered.vcf.gz files in the 
                        parent directory 
  -o OUTPUT, --output OUTPUT
                        directory where to output the filtered or merged files
  -ht, --homozygous_test
                        use argument to collect homozygote statistics
```
