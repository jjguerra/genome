import argparse
import os
import subprocess
import glob


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Format VCF file to use vcftools')
    parser.add_argument('-d', '--directory', '--dir', help='provide parent directory', required=True)
    args = parser.parse_args()

    working_dir = args.directory

    vcf_files = list()
    vcf_files_dir = list()
    vcf_working_dir = list()

    _, family_member_folder_list, _ = os.walk(working_dir).next()

    for member_folder in family_member_folder_list:
        if 'Sample_' in member_folder:
            # complete directory location
            cl_mem_folder = os.path.join(working_dir, member_folder)
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
                    vcf_files.append(member_file)
                    vcf_files_dir.append(member_dir)
                    vcf_working_dir.append(data_dir)

    for vcf_working_dir, vcf_dir, vcf_file in zip(vcf_working_dir, vcf_files_dir, vcf_files):
        print 'processing {0}'.format(vcf_file)

        # remove old filtered file
        for fl in glob.glob(vcf_working_dir + '/*.filtered.vcf.gz'):
            print 'removing {0}'.format(fl)
            os.remove(fl)
        # remove wrong old tbi file
        for fl in glob.glob(vcf_working_dir + '/*vcf.gz.tbi'):
            print 'removing {0}'.format(fl)
            os.remove(fl)

        # create the name for the unzip vcf file
        vcf_filename = vcf_dir.replace('.gz', '')
        subprocess.call(['zcat', vcf_dir, '>', vcf_filename])

        # remove old vcf.gz file
        subprocess.call(['rm -rf', vcf_dir])

        # gz zip the new file with gbzip
        subprocess.call(['bgzip', vcf_filename])

        # tabix the new vcf.gz file
        subprocess.call(['tabix -p vcf', vcf_dir])

        exit(0)