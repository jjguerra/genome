class SNP:

    def __init__(self, base_file, base_snp):

        self.base_snp_file = base_file
        # remove the new line operator
        self.base_snp = base_snp

        # split line by the tab operator
        base_info = base_snp.replace('\n', '').split('\t')

        # make sure values were extracted correctly
        if not base_info:
            raise ValueError('Error while extracting SNPs values on line {0}'.format(base_snp))

        self.chromosome = base_info[0]
        self.position = int(base_info[1])
        self.reference = base_info[2]
        self.one_line_addition = False
        # variable used when adding one line
        self.line = ''
        # variables used when adding more than one line
        self.first_line = ''
        self.second_line = ''

        self.alternate, self.alt_index_dict, self.index_alt_dict = self._allele_index_processing(ref_allele=
                                                                                                 self.reference,
                                                                                                 alt_alleles=
                                                                                                 base_info[3])

        # check if there is more than one genotype information
        if len(base_info[4:]) > 1:
            self.genotype = base_info[4:]
        else:
            self.genotype = list(base_info[-1].replace('\n', ''))

        self.base_written = False
        self.comp_written = False

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
        while len(new_list) != n_col - 1:
            new_list.append('')
        new_list.append(line.split('\t')[-1])
        # return the new line with the tabs operators
        return '\t'.join(new_list)

    @staticmethod
    def _allele_index_processing(ref_allele, alt_alleles):
        """
        This method creates the alternate variable as well as list, dictionaries for its indices and alleles 
        :param ref_allele: variable
        :param alt_alleles: variable or list
        :return: the alternate list, alt_index_dict and index_alt_dict
        """
        # check if there is more than one alternate allele
        if ',' in alt_alleles:
            alternate = alt_alleles.split(',')
        else:
            # make it a list for future analysis
            alternate = list(alt_alleles)

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
        alt_index_dict[ref_allele] = 0
        index_alt_dict[0] = ref_allele

        return alternate, alt_index_dict, index_alt_dict

    @staticmethod
    def _create_line(chrom, pos, ref, alt, gen, n_col):
        """
        Creates the formatted line to be printed
        :param chrom: chromosome
        :param pos: position
        :param ref: reference allele
        :param alt: alternate allele
        :param gen: genotype
        :param n_col: the location of the last genotype
        :return: the created line
        """
        n_alt = ','.join(alt)
        n_gen = '\t'.join(gen[:-1])
        line = '{0}\t{1}\t{2}\t{3}\t{4}'.format(chrom, pos, ref, n_alt, n_gen)
        # need to perform this operation in order to add as tab operands as there are column left
        if len(line) != n_col:
            line += '\t'
        # add the last genotype corresponding to the comp file
        line += gen[-1] + '\n'
        return line

    def add_snp(self, add_snp_file, additional_snp_info, n_col):
        """
        Compares a new SNP information to the previously stored   
        :param add_snp_file:
        :param additional_snp_info: new SNP information 
        """
        # comp_info (list): 0 = chrom, 1 = pos, 2 = ref, 3 = alt, 4 = genotype
        comp_info = additional_snp_info.split('\t')

        # add new SNP information to variables
        n_chrom = comp_info[0]
        n_position = int(comp_info[1])
        n_ref = comp_info[2]
        n_genotype = comp_info[4].replace('\n', '')

        # compare values based on the groups extracted
        # compare chromosome value
        if self.chromosome == n_chrom:

            # if chromosome values are the same, compare position
            if self.position == n_position:

                n_alt, n_alt_index_dict, n_index_alt_dict = self._allele_index_processing(ref_allele=n_ref,
                                                                                          alt_alleles=comp_info[3])

                # make sure REF are the same, if not then error!
                if self.reference != n_ref:
                    msg = 'Error while processing base file {0} and comparison file {1}. \n' \
                          'Base file line = {2} \n' \
                          'Error file line = {3} \n' \
                          'Different REF for the same chromosome and position.'.format(self.base_snp_file,
                                                                                       add_snp_file,
                                                                                       self.base_snp,
                                                                                       additional_snp_info)
                    raise ValueError(msg)

                # list used to keep the genotype data of the base file
                base_genotype_alleles = list()
                # convert base genotype to their alleles
                for base_genotype in self.genotype:
                    if base_genotype:
                        list_base_allele = base_genotype.split('/')
                        first_base_allele = self.index_alt_dict[int(list_base_allele[0])]
                        second_base_allele = self.index_alt_dict[int(list_base_allele[1])]
                        base_genotype_alleles.append(first_base_allele + '/' + second_base_allele)
                    else:
                        base_genotype_alleles.append(base_genotype)

                # split the comp genotype on its alleles
                genotype_indices = n_genotype.split('/')
                # convert comp indices to known alleles
                first_comp_allele = n_index_alt_dict[int(genotype_indices[0])]
                second_com_allele = n_index_alt_dict[int(genotype_indices[1])]
                base_genotype_alleles.append(first_comp_allele + '/' + second_com_allele)

                # add the new ALT allele(s) to the list of base ALT alleles
                self.alternate.extend(n_alt)
                # remove duplicates
                self.alternate = list(set(self.alternate))
                # sort them
                self.alternate.sort()

                # this dictionary has the ALT allele as keys and the ALT allele index as values
                alt_index_dict = dict()
                # this dictionary has the ALT allele index as keys and the ALT allele as values
                index_alt_dict = dict()
                alt_index = 1
                # create the dictionaries based on the types of alternate alleles
                for alt in self.alternate:
                    alt_index_dict[alt] = alt_index
                    index_alt_dict[alt_index] = alt
                    alt_index += 1

                # add the reference information to the dictionaries
                alt_index_dict[self.reference] = 0
                index_alt_dict[0] = self.reference

                # convert all the allele genotype information into indices genotype information
                genotype_list = list()
                for genotype in base_genotype_alleles:
                    if genotype:
                        alleles_list = genotype.split('/')
                        first_allele = alt_index_dict[alleles_list[0]]
                        second_allele = alt_index_dict[alleles_list[1]]
                        allele_line = '{0}/{1}'.format(first_allele, second_allele)
                        genotype_list.append(allele_line)
                    else:
                        genotype_list.append(genotype)

                # create the actual line to go into the file
                self.line = self._create_line(chrom=self.chromosome, pos=self.position, ref=self.reference,
                                              alt=self.alternate, gen=genotype_list, n_col=n_col)

                # change the line flag since only one line will be added to vcf.gz file
                self.one_line_addition = True

                self.base_written = True
                self.comp_written = True

            # if different CHROM and POS values then write on the line whichever is smaller
            elif self.position < n_position:
                self.line = self.base_snp
                self.base_written = True
            else:
                self.line = self._comp_line_format(line=additional_snp_info, n_col=n_col)
                self.comp_written = True

        # if different CHROM then write on the line whichever is smaller
        elif self.chromosome < n_chrom:
            self.line = self.base_snp
            self.base_written = True
        else:
            self.line = self._comp_line_format(line=additional_snp_info, n_col=n_col)
            self.comp_written = True
