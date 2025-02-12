class Sequence:
    """Represents a single sequence in the alignment"""
    def __init__(self, header: str, sequence: str):
        self.header = header
        self.sequence = sequence

    @classmethod
    def from_seqrecord(cls, record) -> 'Sequence':
        """Create a Sequence from a BioPython SeqRecord"""
        return cls(
            header=record.id,
            sequence=str(record.seq)
        )
