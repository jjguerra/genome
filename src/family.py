class Family:
    """
    Family class contains all the structured related to the family members information
    """
    def __init__(self, mother='', father='', son1='', son2='', other=''):
        self.mother = mother
        self.father = father
        self.son1 = son1
        self.son2 = son2
        self.fam5 = other  # this can be an individual value or a list of values

    def get_reference_vcf(self):
        """
        This function return the number of the mother or father to be used as a reference file
        :return: mother file number (if available), father file number (if available) or report an error 
        """
        if self.mother:
            return self.mother
        elif self.father:
            return self.parent
        else:
            raise IOError('No reference file available')
