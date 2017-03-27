import os
import argparse
from vcf import VCF

# current script location
script_dir = '/'.join(os.path.realpath(__file__).split('/')[:-1])


def main():

    # argument method
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', '--dir', help='directory location of the samples')
    parser.add_argument('-f', '--filter', help='use flag to filter columns')
    args = parser.parse_args()

    # check if directory being pass
    if args.directory:
        if not os.path.exists(args.directory):
            raise IOError('directory = "{0}" not found'.format(args.directory))
        else:
            current_directory = args.directory
    else:
        current_directory = script_dir

    # location to store the right directory
    print 'working directory = {0}'.format(current_directory)

    vcf = VCF()
    # read all the vcf files for that family
    vcf.read_files(c_dir=current_directory)

    if args.filter:
        # filter columns of the vcf files
        vcf.filter()

    # merge all the vcf files
    vcf.merge()


if __name__ == '__main__':
    main()
