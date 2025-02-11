class Line:
    def __init__(self, rest_wvl, species, transition, line_width=0.013):
        self.rest_wvl = rest_wvl
        self.species = species
        self.transition = transition
        self.line_width = line_width

        # determine the channel of the lines
        # NIRSpec and MIRI channels and their bounds
        chan_bounds = {
        "nirspec": [2.8708948855637573, 5.269494898093398],
        "ch1-long": [6.530400209798245, 7.649600181524875],
        "ch2-long": [10.010650228883605, 11.699350233480798],
        "ch3-long": [15.41124984738417, 17.978749789996073],
        "ch4-long": [24.40299961855635, 28.69899965589866],
        "ch1-short": [4.900400095357327, 5.739600074157352],
        "ch1-medium": [5.6603998474020045, 6.629999822907848],
        "ch2-short": [7.5106502288836055, 8.770350232312921],
        "ch2-medium": [8.670650076295715, 10.13055008027004],
        "ch3-short": [11.551250190706924, 13.47125014779158],
        "ch3-medium": [13.341250152559951, 15.568750102771446],
        "ch4-short": [17.70300076296553, 20.94900079118088],
        "ch4-medium": [20.693000534083694, 24.47900056699291],
        }

        line_chan = 0
        # loop through the channels to see see if input wvl is within its window
        for chan, bounds in chan_bounds.items():
            if self.rest_wvl > bounds[0] and self.rest_wvl < bounds[1]:
                line_chan = chan
                break  # if the wvl is in the current channel, set line_chan

        self.chan = line_chan

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