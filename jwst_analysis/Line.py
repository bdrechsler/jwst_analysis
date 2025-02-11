class Line:
    def __init__(self, rest_wvl, species, transition, line_width=0.013):
        self.rest_wvl = rest_wvl
        self.sepcies = species
        self.transition = transition
        self.line_width = line_width

    @staticmethod
    def read_line_list(line_list):
        """Functino to read in line list and create line objects

        Args:
            line_list (list of dictionaries): list of lines where each is
            represented by a dictionary. The keys of the dictionary should be
            the same as the attributes of the line class

        Returns:
            return_list (list): list of line objects
        """
        return_list = []
        for item in line_list:
            line = Line(rest_wvl=item['rest_wvl'],
                        species=item['species'],
                        transition=item['transition'],
                        line_width=item['line_width'])
            return_list.append(line)
        return return_list