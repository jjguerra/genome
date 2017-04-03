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
    parser.add_argument('-o', '--output', help='directory where to output the filtered or merged files')
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
    print 'working directory = {0}'.format(current_directory)
    if output_directory:
        print 'output directory = {0}'.format(output_directory)

    vcf = VCF()
    # read all the vcf files for that family
    vcf.read_files(c_dir=current_directory)

    if args.filter:
        # filter columns of the vcf files
        vcf.filter()

    if args.merge:
        # merge all the vcf files
        vcf.merge(output_dir=output_directory)


if __name__ == '__main__':
    main()
