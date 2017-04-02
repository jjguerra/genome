import gzip
from itertools import compress
import os
from family import Family
from shutil import copyfile
from snp import SNP


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

        # these variables are used for merging files
        self.chrom_post_dict = dict()  # key = chrom, val = list of positions

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

        # /Users/jguerra/PycharmProjects/genome/data/Sample_054
        # /Users/jguerra/PycharmProjects/genome/data/Sample_054/Sample_054-
        # extract list of object in the data directory
        # this line only keeps the directories/folder
        _, family_member_file_list, _ = os.walk(self.working_dir).next()

        for member_folder in family_member_file_list:
            # complete directory location
            cl_mem_folder = os.path.join(self.working_dir, member_folder)
            # actual location of the vcf.gz files
            data_dir = os.path.join(cl_mem_folder, 'analysis')
            # get all the files in the analysis location
            # this line only keeps the files
            _, _, member_file_list = os.walk(data_dir).next()

            # check if there exists vcf.gz files or any other files
            # this check that the list is not empty
            if len(member_file_list) == 0:
                raise ValueError('directory = {0} is empty'.format(data_dir))

            for member_file in member_file_list:
                # check for the right vcf.gz file by omitting the vcf.gz.tbi file as well as if a filtered
                # version has already been created
                if member_file.endswith('vcf.gz') and 'filtered' not in member_file:

                    # create a variable with the whole path/directory of the file
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
        Creates a new vcf.gz file with only the CHROM, POS, REF, ALT, genotype columns. This new vcf.gz file will be 
        end in the name 'filtered.vcg.gz'
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

    def _create_base_file(self):
        """
        This function create the base vcf to which to add or merged onto the other vcf files
        :return: the directory of the create vcf reference file and the name 
        """
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

        return ref_file_dir, ref_filename, vcf_filename

    def merge(self, output_dir=''):
        """
        Merge all the filtered.vcf.gz files under the working directory
        """

        print '\nMerge option selected'

        # copy and use the mother's vcf file as a base file
        ref_dir, ref_filename, old_ref_filename = self._create_base_file()

        print '\treference filename = {0}'.format(ref_filename)

        # open the base file
        base_file_obj = gzip.open(ref_dir, 'r')

        # list to keep track of the added files in order to not add duplicates
        merged_files = list()
        # add the reference file
        merged_files.append(old_ref_filename)

        filtered_vcf_files = self.get_vcfs(filtered=True)
        filtered_vcf_files_dirs = self.get_vcfs_dir(filtered=True)

        # loop through all filtered vcf files
        for vcf_dir, vcf_file in zip(filtered_vcf_files_dirs, filtered_vcf_files):

            if vcf_file not in merged_files:

                # open the file to merge
                merging_file_obj = gzip.open(vcf_dir, 'r')
                # open temp merged file
                if output_dir:
                    tmp_file_dir = os.path.join(output_dir, 'tmp.vcf.gz')
                else:
                    tmp_file_dir = ref_dir.replace(ref_filename, 'tmp.vcf.gz')

                tmp_file_obj = gzip.open(tmp_file_dir, 'w+')

                # flag used to skip over the header
                first_iteration = True
                # loop infinitely
                while True:
                    # try and get the next line, if eof then its will just keep looping
                    try:
                        base_line = base_file_obj.next()
                    except StopIteration:
                        pass
                    try:
                        comp_line = merging_file_obj.next()
                    except StopIteration:
                        pass

                    # first iteration contains the header, therefore it needs to be skipped
                    if not first_iteration:

                        # check if not eof for the base file
                        if base_line:
                            # check if not eof for the comp file
                            if comp_line:
                                snp_base = SNP(base_file=ref_filename, base_snp=base_line)
                                snp_base.add_snp(add_snp_file=vcf_file, additional_snp_info=comp_line)

                                # check whether only one line needs to be added to the merged vcf.gz file
                                if snp_base.one_line_addition:
                                    tmp_file_obj.writelines(snp_base.line)
                                # add more than one line
                                else:
                                    tmp_file_obj.writelines(snp_base.first_line)
                                    tmp_file_obj.writelines(snp_base.second_line)

                            else:
                                tmp_file_obj.write(base_line)

                        else:
                            tmp_file_obj.write(comp_line)

                    else:
                        first_iteration = False
                        # use the first line of the reference file
                        line = base_line.replace('\n', '')
                        # obtain the label of the last column of the comparison file
                        new_column = comp_line.split('\t')[-1]
                        # add column name to the reference line
                        line += '\t' + new_column
                        # write line to the temp file
                        tmp_file_obj.writelines(line)

                    # check if eof for both files
                    if not (base_line and comp_line):
                        break

                tmp_file_obj.close()
                merging_file_obj.close()

                # switch the created temp file to the reference file
                copyfile(src=tmp_file_dir, dst=ref_dir)
                # remove tmp file
                os.remove(tmp_file_dir)

                print '\tfinished merging {0}'.format(vcf_file)

        # closing base file
        base_file_obj.close()

        print 'finished merging all vcf files to dir = {0}'.format(ref_filename)
