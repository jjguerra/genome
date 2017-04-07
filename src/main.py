import os
import argparse
from vcf import VCF

# current script location
script_dir = '/'.join(os.path.realpath(__file__).split('/')[:-1])


def main():

    # argument method
    parser = argparse.ArgumentParser()
    # positional argument
    parser.add_argument('-d', '--directory', '--dir', help='parent directory of the samples i.e. Sample_054')
    # optional argument
    parser.add_argument('-f', '--filter', action='store_true',
                        help='use argument to keep the following columns: CHROM, POS, REF, ALT, genotype info')
    # optional argument
    parser.add_argument('-m', '--merge', action='store_true',
                        help='use argument to merge all filtered.vcf.gz files in the parent directory')
    # optional argument
    parser.add_argument('-o', '--output', help='directory where to output the filtered or merged files')
    # optional argument
    parser.add_argument('-ht', '--homozygous_test', action='store_true',
                        help='use argument to collect homozygous statistics')
    # optional argument
    parser.add_argument('-s', '--subset', action='store_true',
                        help='use argument to subset the vcf file based on chromosomes')
    # optional argument
    parser.add_argument('-c', '--chromosome', help='use argument to select the chromosome number on which to subset on')
    # optional argument
    parser.add_argument('-n', '--number_sites', help='use argument to select the number of line on which to subset on',
                        type=int)
    args = parser.parse_args()

    # check if directory being pass
    if args.directory:
        if not os.path.exists(args.directory):
            raise IOError('directory = "{0}" not found'.format(args.directory))
        else:
            current_directory = args.directory
    # instead of making the argument required=True, use this else condition for better
    # error description
    else:
        parser.print_help()
        print ''
        raise IOError('No working directory specified.')

    output_directory = ''
    if args.output:
        if not os.path.exists(args.output):
            raise IOError('output directory = "{0}" not found'.format(args.output))
        else:
            output_directory = args.output

    # location to store the right directory
    print '\nworking directory = {0}'.format(current_directory)
    if output_directory:
        print 'output directory = {0}'.format(output_directory)

    if args.homozygous_test or args.subset:
        hts = True
    else:
        hts = False

    vcf = VCF()
    # read all the vcf files for that family
    vcf.read_files(c_dir=current_directory, homozygous_test_subset=hts)

    if args.filter:
        # filter columns of the vcf files
        vcf.filter()

    if args.merge:
        # merge all the vcf files
        vcf.merge(output_dir=output_directory)

    if args.homozygous_test:
        # collect homozygous statistics
        vcf.homozygous_test(output_dir=output_directory)

    if args.subset:
        if not args.chromosome:
            raise ValueError('Chromosome number must be provided in order to perform subset')

        vcf.subset(chrom=args.chromosome, output_dir=output_directory, n_sites=args.number_sites)

    print 'hheellooo'
if __name__ == '__main__':
    main()
