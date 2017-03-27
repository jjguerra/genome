import gzip
from itertools import compress
import os
from family import Family
from shutil import copyfile


class VCF:
    """
    VCF class contains all the function related to the vcf files i.e. readind, filtering, getting names, etc
    """
    def __init__(self):
        self.working_dir = ''
        self.vcf_family_id = ''
        self.vcf_files = list()
        self.vcf_files_dir = list()
        self.family_info = ''  # this variable will be pointing to the family.Family class

    def read_files(self, c_dir):
        """
        loop through all the folder within the parent directory and stores the vcf.gz files         
        :param c_dir: parent directory 
        """

        # save the c_dir for later use
        self.working_dir = c_dir

        # the number in the Sample_0XX
        self.vcf_family_id = int(c_dir.split('/')[-1].replace('Sample_0', ''))

        # check if c_dir is the right parent directory
        if not self.vcf_family_id:
            raise ValueError('Could not read the family ID from the directory provided')

        # add information to the family_info depending on the c_dir provided
        if self.vcf_family_id == 54:
            self.family_info = Family(mother='054-001', son1='054-003', son2='054-004', other=['054-002',
                                                                                               '054-005'])

        # extract list of object in the data directory
        _, family_member_file_list, _ = os.walk(c_dir).next()

        for member_folder in family_member_file_list:
            # complete directory location
            cl_mem_folder = os.path.join(c_dir, member_folder)
            # actual location of the vcf.gz files
            data_dir = os.path.join(cl_mem_folder, 'analysis')
            # get all the files in the analysis location
            _, _, member_file_list = os.walk(data_dir).next()

            # check if there exists vcf.gz files or any other files
            if len(member_file_list) == 0:
                raise ValueError('Could not get vcf.gz files')

            for member_file in member_file_list:
                # check for the right vcf.gz file by omitting the vcf.gz.tbi file as well as if a filtered
                # version has already been created
                if member_file.endswith('vcf.gz') and 'filtered' not in member_file:

                    member_dir = os.path.join(data_dir, member_file)

                    # add to the object for later processing
                    self.vcf_files.append(member_file)
                    self.vcf_files_dir.append(member_dir)

    def get_vcfs(self, filename='', filtered=''):
        """
        used to get all the vcf
        :param filename: (optional) if filename provided then get that specific vcf filename
        :param filtered: (optional) if true then return name of filtered files
        :return: if filename provided only return one vcf filename else return all the vcf stored
        """

        if filename:
            for vcf_file in self.vcf_files:
                if filename in vcf_file:
                    if filtered:
                        return vcf_file.replace('.vcf.gz', '.filtered.vcf.gz')
                    else:
                        return vcf_file

            # if it did not return a file, then return an error
            raise ValueError('No vcf file with filename={0}'.format(filename))

        # if filtered then modify the names of the files
        if filtered:
            vcf_files_filtered = list()
            for vcf_file in self.vcf_files:
                vcf_files_filtered.append(vcf_file.replace('.vcf.gz', '.filtered.vcf.gz'))

            return vcf_files_filtered
        else:
            return self.vcf_files

    def get_vcfs_dir(self, vcf_file_name='', filtered=''):
        """
        obtain the directory of an specified or all vcf files
        :param vcf_file_name: (optiona) the name of the vcf file to get the interested directory 
        :param filtered: (optional) if true then return name of filtered files
        :return: if vcf_file_name then its directory else all the directories
        """
        # check if asking for a directory of a specific file name
        if vcf_file_name:
            # fetch all the regular or filtered names
            vcf_filenames = self.get_vcfs(filtered=filtered)
            if vcf_file_name not in vcf_filenames:
                raise ValueError('VCF file name "{0}" not found in the list'.format(vcf_file_name))
            else:
                # get the index of interest
                index = vcf_filenames.index(vcf_file_name)
                if filtered:
                    return self.vcf_files_dir[index].replace('.vcf.gz', '.filtered.vcf.gz')
                else:
                    return self.vcf_files_dir[index]

        # if not, then return all the directories
        else:
            if filtered:
                vcf_files_filtered = list()
                for vcf_file in self.vcf_files_dir:
                    vcf_files_filtered.append(vcf_file.replace('.vcf.gz', '.filtered.vcf.gz'))

                return vcf_files_filtered
            else:
                return self.vcf_files_dir

    def filter(self):
        """
        creates a new vcf.gz file with only the CHROM, POS, REF, ALT, genotype columns
        """
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

    def _create_reference_file(self):
        # get the mother or father vcf file number to be used as a reference
        ref_number = self.family_info.get_reference_vcf()
        # get vcf filename
        vcf_filename = self.get_vcfs(filename=ref_number, filtered=True)
        # need to rename the file by removing the reference to the first family member
        ref_filename_list = vcf_filename.split('.')
        ref_filename_list[0] = ref_filename_list[0].split('-')[0]
        ref_filename = '.'.join(ref_filename_list)

        ref_file_dir = os.path.join(self.working_dir, ref_filename)
        # create a reference vcf
        vcf_file_dir = self.get_vcfs_dir(vcf_filename, filtered=True)
        copyfile(vcf_file_dir, ref_file_dir)

        return ref_file_dir, ref_filename

    def merge(self):

        dir, filiname = self._create_reference_file()
