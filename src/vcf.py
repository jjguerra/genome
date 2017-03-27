import os
import argparse

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
                if member_file.endswith('vcf.gz'):

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

if __name__ == '__main__':
    main()
