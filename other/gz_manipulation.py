import argparse
import gzip

"""
Script to write a specified number of lines from a vcf.gz to a new vcf.gz file
"""

# argument method
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--source', help='source vcf.gz file path')
parser.add_argument('-n', '--number_lines', help='number of lines to copy')

args = parser.parse_args()

vcf_dir = args.source

file_obj_read = gzip.open(vcf_dir, 'r')

filename = vcf_dir.split('/')[-1]
new_filename = 'practice_' + filename.split('.')[0] + '.vcf.gz'
new_file_path = vcf_dir.replace(filename, new_filename)

file_obj_write = gzip.open(new_file_path, 'w+')

if args.number_lines:
    limit = int(args.number_lines)
else:
    limit = 10

for index, line in enumerate(file_obj_read):
    import IPython
    IPython.embed()
    print line
    file_obj_write.writelines(line)

    if index == limit:
        break

file_obj_write.close()
file_obj_read.close()

print 'script finished'
print 'file {0} converted to {1}'.format(filename, new_file_path.split('/')[-1])
