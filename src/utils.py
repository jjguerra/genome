import subprocess


def change_file_permission(file_directory=''):
    """
    This function is used in order to change the file permission of any file to o=r
    :param file_directory: 
    """
    subprocess.call(['chmod', '-R', 'o=r', file_directory])


def is_homozygous(genotype):
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


def check_parent_offspring(p_gen, o_gen):
    """
    This function checks the parent and offspring genotype and return the position where the offspring possesses the 
    parent genotype with 0 = False, 1 = True
    :param p_gen: parent genotype
    :param o_gen: offspring genotype
    :return: return a list with two element. First element correspond to the left allele and second element correspond
     to the right allele
    """

    indices_list = [0, 0]

    p_genotypes = p_gen.split('/')
    o_genotypes = o_gen.split('/')

    for gen_indx, gen in enumerate(o_genotypes):
        if gen in p_genotypes:
            indices_list[gen_indx] = 1
        else:
            indices_list[gen_indx] = 0

    if indices_list[0] == 1 and indices_list[1] == 0:
        left_right_both_alleles = 'Left'
    elif indices_list[0] == 0 and indices_list[1] == 1:
        left_right_both_alleles = 'Right'
    elif indices_list[0] == 1 and indices_list[1] == 1:
        left_right_both_alleles = 'Both'
    else:
        left_right_both_alleles = 'None'

    return left_right_both_alleles, indices_list


def index_genotype(ref, alt):
    """
    This method creates the alternate variable as well as list, dictionaries for its indices and alleles 
    :param ref: variable
    :param alt: variable or list
    :return: the alternate list, alt_index_dict and index_alt_dict
    """
    # check if there is more than one alternate allele
    if ',' in alt:
        alternate = alt.split(',')
    else:
        # make it a list for future analysis
        alternate = list(alt)

    # this dictionary has the ALT allele as keys and the ALT allele index as values
    alt_index_dict = dict()
    # this dictionary has the ALT allele index as keys and the ALT allele as values
    index_alt_dict = dict()
    alt_index = 1
    for alt in alternate:
        alt_index_dict[alt] = alt_index
        index_alt_dict[alt_index] = alt
        alt_index += 1

    # add the reference information to the dictionaries
    alt_index_dict[ref] = 0
    index_alt_dict[0] = ref

    return alternate, alt_index_dict, index_alt_dict


def genotype_to_numeric(ref, alt, genotype):
    """
    convert the alphabetic genotype to a numeric genotype
    :param ref: reference alleles
    :param alt: alternative alleles
    :param genotype: alphabetic genotype
    :return: numeric genotype
    """
    _, alt_index_dict, _ = index_genotype(ref, alt)

    new_genotype = list()
    old_genotype = genotype.split('/')
    for og in old_genotype:
        new_genotype.append(alt_index_dict(og))

    return '/'.join(new_genotype)


def genotype_to_alphabetic(ref, alt, genotype):
    """
    convert the numeric genotype to a alphabetic genotype
    :param ref: reference alleles
    :param alt: alternative alleles
    :param genotype: alphabetic genotype
    :return: alphabetic genotype
    """
    _, _, index_alt_dict = index_genotype(ref, alt)

    new_genotype = list()
    old_genotype = genotype.split('/')
    for og in old_genotype:
        new_genotype.append(index_alt_dict(og))

    return '/'.join(new_genotype)

