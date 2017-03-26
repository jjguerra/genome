import os
import argparse

# current script location
script_dir = '/'.join(os.path.realpath(__file__).split('/')[:-1])


def read_files(c_dir):
    """
    loop through all the folder within the parent directory and stores the vcf.gz files 
    
    :param c_dir: parent directory 
    :return: 
    """
    # extract list of object in the data directory
    _, object_list, _ = os.walk(c_dir).next()

    print object_list


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

    read_files(c_dir=current_directory)

    # # extract list of object in the data directory
    # _, object_list, _ = os.walk(data_path).next()

if __name__ == '__main__':
    main()
