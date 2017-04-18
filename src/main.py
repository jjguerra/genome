import os
import argparse
from vcf import VCF

# current script location
script_dir = '/'.join(os.path.realpath(__file__).split('/')[:-1])


def process_arguments(c_directory, arg_parser, o_directory='', chrom='', filter_flag=False, merge_flag=False):
    """
    function takes care of necessary argument information
    :param c_directory: location of the parent directory of the family
    :param arg_parser: use to print the message if error
    :param o_directory: (optional) directory of the output file
    :param chrom: (optional) chromosome being considered
    :param filter_flag: (optional) whether action is to filter the columns of the vcf files
    :param merge_flag: (optional) whether the action is to merge multiple vcf files
    :return: 
    """
    # check if directory being pass
    if c_directory:
        if not os.path.exists(c_directory):
            raise IOError('directory = "{0}" not found'.format(c_directory))
        else:
            working_directory = c_directory
    # instead of making the argument required=True, use this else condition for better
    # error description
    else:
        arg_parser.print_help()
        print ''
        raise IOError('No working directory specified.')

    output_directory = ''
    if o_directory:
        if not os.path.exists(o_directory):
            raise IOError('output directory = "{0}" not found'.format(o_directory))
        else:
            output_directory = o_directory

    # location to store the right directory
    print '\nworking directory = {0}'.format(working_directory)
    if output_directory:
        print 'output directory = {0}'.format(output_directory)

    if chrom:
        chromosome = chrom
    else:
        chromosome = ''

    # check if a filtering or merging is going to be performed
    filter_merge_flag = False
    if filter_flag or merge_flag:
        filter_merge_flag = True

    return working_directory, output_directory, chromosome, filter_merge_flag


def main():

    # argument method
    parser = argparse.ArgumentParser()
    # positional argument
    parser.add_argument('-d', '--directory', '--dir', help='parent directory of the samples i.e. Sample_054',
                        required=True)
    # optional argument
    parser.add_argument('-f', '--filter', action='store_true',
                        help='use argument to keep the following columns: CHROM, POS, REF, ALT, genotype info. \n'
                             'Note: if option selection, the -d can either be the vcf.gz or the parent directory. '
                             'If vcf.gz, it will perform filtering on that file. Else, it will perform filtering on'
                             'all the files in the hierarchical directory.')
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
    # optional argument
    parser.add_argument('-p', '--phase', action='store_true', help='use argument to select the phasing test')
    args = parser.parse_args()

    working_directory, output_directory, chromosome, filter_merge = process_arguments(c_directory=args.directory,
                                                                                      o_directory=args.output,
                                                                                      filter_flag=args.filter,
                                                                                      merge_flag=args.merge,
                                                                                      chrom=args.chromosome,
                                                                                      arg_parser=parser)

    # create VCF object
    vcf = VCF()

    # filtering and merging have to read all the vcf files on the subdirectories of the families
    if filter_merge:
        # read all the vcf files for that family
        vcf.read_files(c_dir=working_directory)
        # filter columns of the vcf files
        if args.filter:
            vcf.filter(output_dir=output_directory)
        else:
            vcf.merge(output_dir=output_directory)

    # subset, homozygous or phasing tests have to read a vcf file in the parent directory
    else:
        # for the subset, the chromosome should not be given in the read_file function
        if args.subset:
            # chromosome if required if goal is to subset since the file will be subsetted on it
            if not args.chromosome:
                raise ValueError('Chromosome number must be provided in order to perform subset')

            # read all the vcf files for that family
            vcf.read_files(c_dir=working_directory, homozygous_test_subset=True)
            vcf.subset(chrom=args.chromosome, output_dir=output_directory, n_sites=args.number_sites)

        else:
            # read all the vcf files for that family
            vcf.read_files(c_dir=working_directory, homozygous_test_subset=True, chrom=chromosome)

            if args.homozygous_test:
                # collect homozygous statistics
                vcf.tests(output_dir=output_directory, chrom=chromosome, homozygous_test=True)

            if args.phase:
                vcf.tests(output_dir=output_directory, chrom=chromosome, homozygous_test=False)


if __name__ == '__main__':
    main()
