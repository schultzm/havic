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

    def __init__(self, alignment, refisolate, gap_char):

        self.alignment = alignment
        self.refisolate = refisolate
        self.gap_char = gap_char

    def _get_isolate_coords(self):
        """
        Get the coords of the sequence excluding the 5' and 3' gap padding.

        """
        start_pos = None
        end_pos = None
        for seq in self.alignment:
            if seq.id == self.refisolate:
                print(str(seq.seq))

                from re import finditer
                for match in finditer(f"{self.gap_char}[A-Z]",
                                      str(seq.seq).upper()):
                    start_pos = (match.span())
                for match in finditer(f"{self.gap_char}[A-Z]",
                                      str(seq.seq).upper()[::-1]):
                    end_pos = (-match.end(), -match.start())
                    # end_pos = (match.span(), match.group())
                print(seq.seq[start_pos[0] + 1:end_pos[0] + 1])
                break

                # print('reverse', start_pos, end_pos)
            # break

        # site_set = {"-"}
        # start_pos = 0
        # while len(site_set) == 1:
        #     site_set = set(alignment[:, start_pos])
        #     start_pos += 1
        # site_set = {"-"}
        # # subtract one due to 0 and 1-based indexing issues
        # end_pos = alignment.get_alignment_length() - 1
        # while len(site_set) == 1:
        #     site_set = set(alignment[:, end_pos])
        #     end_pos -= 1
        #
        #
        # pass

        # print(dir(self.alignment))
