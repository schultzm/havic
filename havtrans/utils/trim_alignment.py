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

        We'll use the following MSA as an example:
        >>> from Bio.Alphabet import generic_dna
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("---AAAACGT--A", generic_dna), id="Alpha")
        >>> b = SeqRecord(Seq("---AAA-CGT---", generic_dna), id="Beta")
        >>> c = SeqRecord(Seq("---AAAAGGT---", generic_dna), id="Gamma")
        >>> d = SeqRecord(Seq("---AAAACGT---", generic_dna), id="Delta")
        >>> e = SeqRecord(Seq("---AAA-GGT---", generic_dna), id="Epsilon")
        >>> align = MultipleSeqAlignment([a, b, c, d, e], generic_dna)
        >>> print(align)
        DNAAlphabet() alignment with 5 rows and 13 columns
        ---AAAACGT--A Alpha
        ---AAA-CGT--- Beta
        ---AAAAGGT--- Gamma
        ---AAAACGT--- Delta
        ---AAA-GGT--- Epsilon
        >>> refiso = Trimmed_alignment(align, 'Beta').isolate
        >>> print(refiso)
        Beta
        >>> for record in align:
        ...     if record.id == refiso:
        ...         print(record)
        ID: Beta
        Name: <unknown name>
        Description: <unknown description>
        Number of features: 0
        Seq('---AAA-CGT---', DNAAlphabet())
        

        """
        pass

        # print(dir(self.alignment))
