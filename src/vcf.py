import os
import argparse
import gzip
from itertools import compress


# current script location
script_dir = '/'.join(os.path.realpath(__file__).split('/')[:-1])


class VCF:
    def __init__(self):
        self.vcf_files = list()
        self.vcf_files_dir = list()

    def read_files(self, c_dir):
        """
        loop through all the folder within the parent directory and stores the vcf.gz files         
        :param c_dir: parent directory 
        """
        # extract list of object in the data directory
        _, family_member_file_list, _ = os.walk(c_dir).next()

        for member_folder in family_member_file_list:
            # complete directory location
            cl_mem_folder = os.path.join(c_dir, member_folder)
            # actual location of the vcf.gz files
            data_dir = os.path.join(cl_mem_folder, 'analysis')
            # get all the files in the analysis location
            _, _, member_file_list = os.walk(data_dir).next()
            for member_file in member_file_list:
                # check for the right vcf.gz file by omitting the vcf.gz.tbi file as well as if a filtered
                # version has already been created
                if member_file.endswith('vcf.gz') and 'filtered' not in member_file:

                    member_dir = os.path.join(data_dir, member_file)

                    # add to the object for later processing
                    self.vcf_files.append(member_file)
                    self.vcf_files_dir.append(member_dir)

    def get_vcfs(self):
        return self.vcf_files

    def get_vcfs_dir(self, vcf_file_name=''):
        """
        obtain the directory of an specified or all vcf files
        :param vcf_file_name: (optiona) the name of the vcf file to get the interested directory 
        :return: if vcf_file_name then its directory else all the directories
        """
        # check if asking for a directory of a specific file name
        if vcf_file_name:
            if vcf_file_name not in self.vcf_files:
                raise ValueError('VCF file name "{0}"not found in the list'.format(vcf_file_name))
            else:
                # get the index of interest
                index = self.vcf_files.index(vcf_file_name)
                return self.vcf_files_dir[index]

        # if not, then return all the directories
        else:
            return self.vcf_files_dir

    def filter(self):

        for vcf_dir, vcf in zip(self.vcf_files_dir, self.vcf_files):
            print 'processing {0}'.format(vcf)

            # this flag indicates when to start collecting data
            val_flag = False
            # only take = CHROM, POS, REF, ALT, genotype
            values_of_interest = [True, True, False, True, True, False, False, False, False, True]

            file_obj_read = gzip.open(vcf_dir, 'r')
            # new file filtered
            new_directory = vcf_dir.replace('.vcf.gz', '.filtered.vcf.gz')
            file_obj_write = gzip.open(new_directory, 'w+')
            # genotype title has to be different depending on the vcf file
            genotype = vcf.split('.')[0]
            # file header
            file_obj_write.writelines('#CHROM\tPOS\tREF\tALT\t' + genotype + '\n')

            header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
            for line in file_obj_read:
                if header in line:
                    val_flag = True
                elif val_flag:
                    split_line = line.split('\t')
                    # filter the list split_line based on the boolean list values_of_intest
                    new_list = list(compress(split_line, values_of_interest))
                    # only get the genotype information - remove the metadata after the genotype information
                    new_list[4] = new_list[4].split(':')[0]
                    new_line = '\t'.join(new_list)
                    new_line += '\n'
                    file_obj_write.writelines(new_line)

            # closing files
            file_obj_read.close()
            file_obj_write.close()


def main():

    # argument method
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', '--dir', help='directory location of the samples')
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
    # filter the right column for the vcf files
    vcf.filter()

if __name__ == '__main__':
    main()
