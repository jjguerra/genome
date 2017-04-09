import gzip
from itertools import compress
import os
from family import Family
from shutil import copyfile
from snp import SNP
import subprocess


class VCF:
    """
    VCF class contains all the function related to the vcf files i.e. reading, filtering, getting names, etc
    """

    def __init__(self):
        self.working_dir = ''
        self.vcf_family_id = ''
        self.vcf_files = list()
        self.vcf_files_dir = list()
        self.family_info = Family()  # this variable will be pointing to the family.Family class

    def read_files(self, c_dir, homozygous_test_subset=False, chrom=''):
        """
        loop through all the folder within the parent directory and stores the vcf.gz files         
        :param c_dir: parent directory 
        :param homozygous_test_subset: flag in order to collect a different kind of data for statistics collection
        :param chrom: chromosome number if need of specific file with that chromosome
        """
        # save the c_dir for later use
        self.working_dir = c_dir

        # the number in the Sample_0XX
        self.vcf_family_id = int(c_dir.split('/')[-1].replace('Sample_', ''))

        # check if c_dir is the right parent directory
        if not self.vcf_family_id:
            raise ValueError('Could not read the family ID from the directory provided')

        # add information to the family_info depending on the c_dir provided
        if self.vcf_family_id == 54:
            self.family_info = Family(parent='054-001', offspring=['054-003', '054-004'], other=['054-002', '054-005'])
        elif self.vcf_family_id == 85:
            self.family_info = Family(parent='085-001', offspring=['085-002', '085-006'], other=['085-004', '054-005'])
        elif self.vcf_family_id == 89:
            self.family_info = Family(parent='089-001', offspring=['089-007'], other=['089-005', '089-006', '089-009'])
        elif self.vcf_family_id == 95:
            pass
        elif self.vcf_family_id == 97:
            pass
        elif self.vcf_family_id == 109:
            self.family_info = Family(parent='109-003', offspring=['109-001', '109-002'], other=['109-004'])
        elif self.vcf_family_id == 110:
            pass
        elif self.vcf_family_id == 115:
            self.family_info = Family(parent='115-004', offspring=['115-001', '115-002', '115-003',
                                                                   '115-005', '115-006'])

        # check if chromosome number was passed
        if chrom:
            # add the right definition for file name
            chrom = 'chrom' + str(chrom)

        if not homozygous_test_subset:
            # /Users/jguerra/PycharmProjects/genome/data/Sample_054
            # /Users/jguerra/PycharmProjects/genome/data/Sample_054/Sample_054-
            # extract list of object in the data directory
            # this line only keeps the directories/folder
            _, family_member_folder_list, _ = os.walk(self.working_dir).next()

            for member_folder in family_member_folder_list:
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

        # for stats collection, only look for the cvf.gz file
        else:
            # get all the files in the given working directory
            _, _, dir_file_list = os.walk(self.working_dir).next()

            # loop through all the files
            for vcf_file in dir_file_list:
                # look for the filtered.vcf.gz with the specified family id
                if ('filtered' in vcf_file) and (str(self.vcf_family_id) in vcf_file):
                    # check if the chromosome was pass
                    if chrom:
                        # if chromosome was provided, make sure that the file has it
                        if chrom in vcf_file:
                            # add its information to the class variables
                            self.vcf_files = vcf_file
                            self.vcf_files_dir = os.path.join(self.working_dir, vcf_file)

                    else:
                        # add its information to the class variables
                        self.vcf_files = vcf_file
                        self.vcf_files_dir = os.path.join(self.working_dir, vcf_file)

    def get_vcfs(self, filename='', filtered=False):
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

    def get_vcfs_dir(self, vcf_file_name='', filtered=False):
        """
        obtain the directory of an specified or all vcf files
        :param vcf_file_name: (optional) the name of the vcf file to get the interested directory 
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
        print '\nFilter option selected'
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

            self._change_file_permission(file_directory=new_directory)

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

        self._change_file_permission(file_directory=ref_file_dir)

        return ref_file_dir, ref_filename, vcf_filename

    @staticmethod
    def _comp_line_format(line, n_col):
        """
        This function create the right string to be printed if the comp line is used i.e. it prints the genotype
        information of the comp line in the right column
        :param line: line to be printed
        :param n_col: column index of the genotype index
        :return: the new line with the right spaces
        """
        new_list = line.split('\t')[:4]
        while len(new_list) != n_col:
            new_list.append('')
        new_list.append(line.split('\t')[-1])
        # return the new line with the tabs operators
        return '\t'.join(new_list)

    def merge(self, output_dir=''):
        """
        Merge all the filtered.vcf.gz files under the working directory
        """

        print '\nMerge option selected'

        # copy and use the parent's vcf file as a base file
        ref_dir, ref_filename, old_ref_filename = self._create_base_file()

        print '\treference filename = {0}'.format(ref_filename)

        # list to keep track of the added files in order to avoid duplicates
        merged_files = list()
        # add the reference file
        merged_files.append(old_ref_filename)

        # get all the vcf and their location with the filtered flag
        filtered_vcf_files = self.get_vcfs(filtered=True)
        filtered_vcf_files_dirs = self.get_vcfs_dir(filtered=True)

        n_file = 0
        # loop through all filtered vcf files
        for vcf_dir, vcf_file in zip(filtered_vcf_files_dirs, filtered_vcf_files):

            n_file += 1

            # check the vcf is not in the merged files
            if vcf_file not in merged_files:

                # open the base file
                base_file_obj = gzip.open(ref_dir, 'r')
                # open the file to merge
                merging_file_obj = gzip.open(vcf_dir, 'r')

                # create a temporary merged file
                if output_dir:
                    tmp_file_dir = os.path.join(output_dir, 'tmp.vcf.gz')
                else:
                    tmp_file_dir = ref_dir.replace(ref_filename, 'tmp.vcf.gz')

                # open temp merged file
                tmp_file_obj = gzip.open(tmp_file_dir, 'w+')

                # flag used to skip over the header
                header = True

                # initialize lines variables
                base_line = ''
                comp_line = ''

                # This number will be used to know where to print comp information i.e. information used to print the
                # comp column in the right order
                comp_column = 0

                # loop infinitely
                while True:

                    # flags used to turn on deletion of their respective lines based on whether the line was
                    # written to the temp file or not
                    reset_base_line = False
                    reset_comp_line = False

                    # only read the file when the variable is empty
                    if not base_line:
                        # try and get the next line, if eof then it will just keep looping
                        try:
                            base_line = base_file_obj.next()
                        except StopIteration:
                            pass
                    # only read the file when the variable is empty
                    if not comp_line:
                        # try and get the next line, if eof then it will just keep looping
                        try:
                            comp_line = merging_file_obj.next()
                        except StopIteration:
                            pass

                    # check if reading the header line
                    if header:
                        # turn of the header flag
                        header = False
                        # use the first line of the reference file and remove its new line operator
                        line = base_line.replace('\n', '')
                        # obtain the name of the last column of the comparison file i.e. the name of the person
                        new_column = comp_line.split('\t')[-1]
                        # add column name to the reference line
                        line += '\t' + new_column
                        # write line to the temp file
                        tmp_file_obj.writelines(line)
                        # obtain the number of column in the reference line. This number will be used to know where
                        # to print comp information
                        comp_column = len(line.split('\t'))

                        # flags used to delete the content of their respective line after processing
                        reset_base_line = True
                        reset_comp_line = True

                    else:
                        # check if not eof for the base file
                        if base_line:
                            # check if not eof for the comp file
                            if comp_line:
                                snp_base = SNP(base_file=ref_filename, base_snp=base_line)
                                snp_base.add_snp(add_snp_file=vcf_file, additional_snp_info=comp_line,
                                                 n_col=comp_column)

                                # always add a line based on the following requirements
                                #   - chrom and position are the same
                                #   - the information with the lowest chrom and position
                                tmp_file_obj.writelines(snp_base.line)

                                # turn on the respective flags
                                if snp_base.base_written:
                                    reset_base_line = True
                                if snp_base.comp_written:
                                    reset_comp_line = True

                            else:
                                # write the base information to the tmp file
                                tmp_file_obj.write(base_line)
                                # delete the content in the base line
                                reset_base_line = True

                        else:
                            # check if not eof for the comp file
                            if comp_line:
                                # use the right format to print the line i.e. place the genotype in the right column
                                new_line = self._comp_line_format(line=comp_line, n_col=comp_column)
                                # write the new line to the temp file
                                tmp_file_obj.writelines(new_line)
                                # delete the content in the comp line
                                reset_comp_line = True

                    # check if eof for both files
                    if not base_line and not comp_line:
                        break

                    # if these flags have been turned on by having written those lines in the file, then remove
                    # the values in the lines.
                    if reset_base_line:
                        base_line = ''
                    if reset_comp_line:
                        comp_line = ''

                # close all the files
                tmp_file_obj.close()
                merging_file_obj.close()
                base_file_obj.close()

                # switch the created temp file to the reference file
                copyfile(src=tmp_file_dir, dst=ref_dir)
                # remove tmp file
                os.remove(tmp_file_dir)

                print '\tfinished merging {0}'.format(vcf_file)

        print 'finished merging all vcf files to dir = {0}'.format(ref_filename)

    @staticmethod
    def _is_homozygous(genotype):
        """
        check if the genotype passed is homozygous
        :param genotype: genotype info i.e. 0/0 
        :return: true if it is homozygous else false
        """
        genotype_alleles = genotype.split('/')

        if genotype_alleles[0] == genotype_alleles[1]:
            return True
        else:
            return False

    def homozygous_test(self, output_dir='', chrom=''):
        """
        This function collect homozygous statistics from mother to offsprings
        :param output_dir: (optional) location to output the statistics file
        :param chrom: chromosome being evaluated
        """
        print '\nhomozygous_test option selected'

        # open the vcf file
        file_obj_read = gzip.open(self.vcf_files_dir, 'r')

        # create the statistics filename
        if chrom:
            filename = 'homozygous_test_' + str(self.vcf_family_id) + '_' + 'chrom' + chrom + '.text'
        else:
            filename = 'homozygous_test_' + str(self.vcf_family_id) + '.text'

        # check if output directory was provided
        if output_dir:
            homozygous_dir = os.path.join(output_dir, filename)
        # if not provided, use the working directory
        else:
            homozygous_dir = os.path.join(self.working_dir, filename)

        # output object
        file_obj_write = open(homozygous_dir, 'w+')

        # keep track of different alleles between homozygous parent and offspring
        mismatch_parent_offspring_homo = [0.0] * len(self.family_info.offspring)
        # keep track of different alleles between parent and offspring
        mismatch_parent_offspring_all = [0.0] * len(self.family_info.offspring)
        # total numbers of sites
        total_num_sites = float(0)
        # total number of sites evaluated
        total_num_sites_eval = float(0)

        # set flag to skip the header
        header = True

        # keep the parent column index
        parent_col_index = ''

        # read the file line by line
        for line in file_obj_read:

            # skip and write the header
            if header:
                header = False
                header_columns = line.split('\t')
                file_obj_write.writelines(line)

                parent_number = self.family_info.parent

                # loop through all the element of the header_columns to find the right col index of the parent
                for header_index, header_info in enumerate(header_columns):
                    if parent_number in header_info:
                        parent_col_index = header_index
                        break

                if parent_col_index == '':
                    raise ValueError('The parent was not found in the vcf file')

                # obtain a list of offspring IDs
                offspring_list = self.family_info.offspring

                # for each offspring, find their column in the vcf.gz file and make a dictionary where the key
                # is the offspring and its value is the column
                offspring_col_index_dict = dict()
                for offspring in offspring_list:
                    for index, columns in enumerate(header_columns):
                        if offspring in columns:
                            offspring_col_index_dict[offspring] = index
                            break

                if len(offspring_col_index_dict.keys()) != len(offspring_list):
                    # Note: we can actually check for the offsprings which values were not populated
                    raise ValueError('Some of the offspring were not found in the vcf file')

            # starting from the second line
            else:
                total_num_sites += 1
                # parse the line
                line_information = line.split('\t')

                # parse the genotypes of parent
                genotype_parent = line_information[parent_col_index]
                # use for a test later on
                genotype_parent_list = line_information[parent_col_index].split('/')

                # check if the parent's genotype was provided since there are cases when it is not.
                # if it is then check whether it is homozygous
                if genotype_parent and self._is_homozygous(genotype=genotype_parent):
                    # increased sites evaluated
                    total_num_sites_eval += 1

                    for offspring_index, offspring in enumerate(offspring_list):
                        try:
                            # get the offspring column index
                            col_index = offspring_col_index_dict[offspring]

                            # get the genotype information
                            offspring_genotype = line_information[col_index].replace('\n', '')

                            # check if offspring is homozygous
                            if self._is_homozygous(offspring_genotype):
                                # if the parent and the child are homozygous on different allele,
                                # print a warning message
                                # write the line into the file
                                if offspring_genotype[0] != genotype_parent[0]:
                                    mismatch_parent_offspring_homo[offspring_index] += 1
                                    file_obj_write.writelines(line)
                                    print 'error on chrom = {0}, position = {1} - parents and offspring homozygous ' \
                                          'alleles mismatch'.format(line_information[0], line_information[1])
                                    print 'offspring = {0}'.format(offspring)
                                    print line
                            # if its not homozygous, make sure it has at least one allele transmitted by the parent
                            else:
                                first_allele, second_allele = offspring_genotype.split('/')
                                if (first_allele not in genotype_parent_list) and\
                                        (second_allele not in genotype_parent_list):
                                    mismatch_parent_offspring_all[offspring_index] += 1
                                    file_obj_write.writelines(line)
                                    msg = 'error chrom = {0}, position = {1} - parents and offspring alleles mismatch'.\
                                        format(line_information[0], line_information[1])
                                    print msg
                                    print 'offspring = {0}'.format(offspring)
                                    print line

                        # this is done for sites where the offspring does not have information
                        except IndexError:
                            pass

        msg = '\ntotal number of sites = {0}'.format(total_num_sites)
        self._output_line(file_obj=file_obj_write, line_info=msg)
        msg = 'total number of sites evaluated = {0}'.format(total_num_sites_eval)
        self._output_line(file_obj=file_obj_write, line_info=msg)
        # make a line
        self._output_line(file_obj=file_obj_write, line_info='')

        for offspring_index, offspring in enumerate(self.family_info.offspring):

            # provide offspring information
            msg = 'offspring = {0}'.format(offspring)
            self._output_line(file_obj=file_obj_write, line_info=msg)

            # provide information on homozygous parent and offspring difference
            msg = 'number of wrong homozygous parent and homozygous offspring SNP = {0}'.format(
                mismatch_parent_offspring_homo[offspring_index])
            self._output_line(file_obj=file_obj_write, line_info=msg)

            ratio = mismatch_parent_offspring_homo[offspring_index]/total_num_sites_eval
            msg = 'mismatch ratio = {0}'.format(ratio)
            self._output_line(file_obj=file_obj_write, line_info=msg)

            # provide information about homozygous parent and different offspring alleles
            msg = 'number of wrong alleles between parent and offspring SNP = {0}'.format(
                mismatch_parent_offspring_all[offspring_index])
            self._output_line(file_obj=file_obj_write, line_info=msg)

            ratio = mismatch_parent_offspring_all[offspring_index]/total_num_sites_eval
            msg = 'mismatch ratio = {0}'.format(ratio)
            self._output_line(file_obj=file_obj_write, line_info=msg)
            self._output_line(file_obj=file_obj_write, line_info='\n')

        print '{0} written\n'.format(self.vcf_files)

        # close the files and print messages
        file_obj_write.close()
        file_obj_read.close()

        self._change_file_permission(file_directory=homozygous_dir)

    @staticmethod
    def _output_line(file_obj, line_info):
        print line_info
        file_obj.writelines(line_info + '\n')

    def subset(self, chrom, output_dir='', n_sites=float("inf")):
        """
        This function subsets the vcf.gz file based on the chromosome number
        :param chrom: chromosome
        :param output_dir: (optional) location to output the statistics file
        :param n_sites: the total number of sites to include in the file
        """

        print '\nsubset option selected on chrom {0}'.format(chrom)

        # open the vcf file
        file_obj_read = gzip.open(self.vcf_files_dir, 'r')

        # create the output filename
        subset_filename = 'chrom{0}.vcf.gz'.format(chrom)
        filename = self.vcf_files_dir.replace('vcf.gz', subset_filename)

        # check if output directory was provided
        if output_dir:
            subset_dir = os.path.join(output_dir, filename)
        # if not provided, use the working directory
        else:
            subset_dir = os.path.join(self.working_dir, filename)

        # output object
        file_obj_write = gzip.open(subset_dir, 'w+')

        # keep track of sites added
        num_sites_added = 0

        header = True
        for line in file_obj_read:

            if header:
                header = False
                file_obj_write.writelines(line)

            else:
                # split the lab on the tab operator and grab the chromosome number
                line_chromosome = line.split('\t')[0]
                if line_chromosome == str(chrom):
                    file_obj_write.writelines(line)
                    num_sites_added += 1
                    if num_sites_added == n_sites:
                        break

        print 'total number of sites included = {0}'.format(num_sites_added)
        print 'finished obtaining subset of the file\n'

        file_obj_read.close()
        file_obj_write.close()

        self._change_file_permission(file_directory=subset_dir)

    @staticmethod
    def _change_file_permission(file_directory=''):
        """
        This function is used in order to change the file permission of any file to o=r
        :param file_directory: 
        """
        subprocess.call(['chmod', '-R', 'o=r', file_directory])
