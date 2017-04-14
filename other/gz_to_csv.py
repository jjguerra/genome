import argparse
import gzip
import os


def main():

    # argument method
    parser = argparse.ArgumentParser()
    # positional argument
    parser.add_argument('-d', '--directory', '--dir', help='parent directory of the samples', required=True)
    args = parser.parse_args()

    working_dir = args.directory
    _, _, files_under_wd = os.walk(working_dir).next()

    for vcf_file in files_under_wd:
        if 'filtered.vcf.gz' in vcf_file:
            working_dir_file = os.path.join(working_dir, vcf_file)
            break

    file_obj_read = gzip.open(working_dir_file, 'r')

    csv_filename = working_dir_file.replace('vcf.gz', 'vcf.csv')
    file_obj_write = open(csv_filename, 'w+')

    for line in file_obj_read:
        n_line = line.replace(',', '|')
        n_line = n_line.replace('\t', ',')
        file_obj_write.write(n_line)

if __name__ == '__main__':
    main()

