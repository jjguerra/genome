class SNP:
    def __init__(self, base_file, base_snp):

        self.base_snp_file = base_file
        # remove the new line operator
        self.base_snp = base_snp.replace('\n', '')

        # split line by the tab operator
        base_info = base_snp.split('\t')

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

        # check if there is more than one alternative allele
        if ',' in base_info[3]:
            self.alternate = base_info[3].split(',')
        else:
            # make it a list for future analysis
            self.alternate = [base_info[3]]

        # check if there is more than one genotype information
        if len(base_info[4:]) > 1:
            self.genotype = base_info[-1].replace('\n', '')
        else:
            self.genotype = [base_info[-1].replace('\n', '')]

    def add_snp(self, add_snp_file, additional_snp_info):
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
        n_alt = comp_info[3]
        n_gen = comp_info[4]

        # compare values based on the groups extracted
        # compare chromosome value
        if self.chromosome == n_chrom:

            # if chromosome values are the same, compare position
            if self.position == n_position:

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

                # compare ALT alleles
                if n_alt in self.alternate:
                    # NOTE --- CHECK THAT ALT alleles MATCH genotype info
                    # if same ALT alleles, add the new genotype information to line
                    self.line = self.base_snp + n_gen + '\n'
                    # change the line flag since only one line will be added to vcf.gz file
                    self.one_line_addition = True

                else:
                    # if alleles are different
                    #   1. need to add the new ALT allele to the ALT column in the right order
                    #   2. need to modify the right genotype to account for the new allele
                    # add diff allele

                    # allele_index in the possible genotype number for the allele
                    # REF allele is 0
                    allele_dict = dict()
                    allele_dict[self.reference] = 0

                    allele_index = 1
                    for allele in self.alternate:
                        allele_dict[allele] = allele_index
                        allele_index += 1

                    # add the new ATL allele to the dict of ALT alleles with its respective index
                    allele_dict[n_alt] = allele_index
                    # add the new ATL allele to the list of ALT alleles
                    self.alternate.append(n_alt)
                    # sort the ALT allele list since it needs to be sorted to be added to the ALT column
                    self.alternate = self.alternate.sort()

                    for c_genotype in self.genotype:
                        pass

            # if different POS values then whichever is smaller, write it first
            elif self.position < n_position:
                self.first_line = self.base_snp
                self.second_line = additional_snp_info
            else:
                self.first_line = additional_snp_info
                self.second_line = self.base_snp

        # if different chromosome numbers then whichever is smaller, write it first
        elif self.chromosome < n_chrom:
            self.first_line = self.base_snp
            self.second_line = additional_snp_info
        else:
            self.first_line = additional_snp_info
            self.second_line = self.base_snp
