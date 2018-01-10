"""A class for trimming a BioPython MultipleSeqAlignment object.

Inherits from BioPython MultipleSeqAlignment

Input:
    MultipleSeqAlignment
Output:
    Trimmed MultipleSeqAlignment.
"""

from Bio.Align import MultipleSeqAlignment


class Trimmed_alignment():
    """Trim a BioPython MSA object.

    Given an alignment, trim the 5' and 3' gap-only regions.
    """

    def __init__(self, alignment, refisolate):
        self.alignment = alignment
        self.isolate = refisolate

    def _get_isolate_coords(self):
        """
        Get the coords of the sequence excluding the 5' and 3' gap padding.
        """
        print(self.alignment[0])


        # print(dir(self.alignment))
