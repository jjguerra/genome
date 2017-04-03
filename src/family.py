class Family:
    """
    Family class contains all the structured related to the family members information
    """
    def __init__(self, parent='', offspring='', other=''):
        self.parent = parent
        self.offspring = offspring
        self.fam5 = other  # this can be an individual value or a list of values

    def get_reference_vcf(self):
        """
        This function return the id of the parent to be used as a reference file
        :return: parent file number (if available), father file number (if available) or report an error 
        """
        # check if parent was recorded
        if self.parent:
            # return parents id
            return self.parent
        else:
            raise IOError('No reference file available')
