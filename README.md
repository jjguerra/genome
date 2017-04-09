# genome
```
usage: main.py [-h] [-d DIRECTORY] [-f] [-m] [-o OUTPUT] [-ht] [-s]
               [-c CHROMOSOME] [-n NUMBER_SITES]

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
                        use argument to collect homozygous statistics
  -s, --subset          use argument to subset the vcf file based on
                        chromosomes
  -c CHROMOSOME, --chromosome CHROMOSOME
                        use argument to select the chromosome number on which
                        to subset on
  -n NUMBER_SITES, --number_sites NUMBER_SITES
                        use argument to select the number of line on which to
                        subset on
```

Use the following commands depending on the function.

Filter:
```
python main.py -f -d DIRECTORY [-o OUTPUT DIRECTORY]
```
Merge:
```
python main.py -m -d DIRECTORY [-o OUTPUT DIRECTORY]
```
Subset:
```
python main.py -s -c CHROMOSOME -d DIRECTORY [-o OUTPUT DIRECTORY] [-n NUMBER SITES per CHROMOSOME]
```
Homozygous Test:
```
python main.py -ht [-c CHROMOSOME] -d DIRECTORY [-o OUTPUT DIRECTORY]
```




