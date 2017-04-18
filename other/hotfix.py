import gzip
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Format VCF file to use vcftools')
    parser.add_argument('-d', '--directory', '--dir', help='provide parent directory', required=True)
    args = parser.parse_args()

    working_dir = args.directory

    file_obj_read = gzip.open(working_dir, 'r')

    n_filename = working_dir.replace('.vcf.gz', 'fixed.vcf.gz')
    file_obj_write = gzip.open(n_filename, 'w+')

    for line in file_obj_read:
        if '#CHROM\tPOS\tID' in line:
            import IPython
            IPython.embed()
            file_obj_write.writelines(line)
        else:
            file_obj_write.writelines(line)

    file_obj_read.close()
    file_obj_write.close()
