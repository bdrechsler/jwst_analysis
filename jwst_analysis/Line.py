class Line:
    def __init__(self, rest_wvl, species, transition, line_width=0.013):
        self.rest_wvl = rest_wvl
        self.sepcies = species
        self.transition = transition
        self.line_width = line_width