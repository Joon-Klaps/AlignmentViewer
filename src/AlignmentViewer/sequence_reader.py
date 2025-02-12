from typing import List, Union, TextIO
from pathlib import Path
from Bio import SeqIO
from .sequence import Sequence

class SequenceReader:
    """Handles sequence file parsing using BioPython"""

    @staticmethod
    def parse(source: Union[str, Path, TextIO], max_seqs: int) -> List[Sequence]:
        """Parse sequence format from various input sources using BioPython SeqIO"""
        try:
            if isinstance(source, (str, Path)):
                format = SequenceReader._guess_format(str(source))
                with open(source) as f:
                    records = list(SeqIO.parse(f, format))
            else:
                records = list(SeqIO.parse(source, "fasta"))

            max_seqs = 0 if max_seqs > len(records) else max_seqs
            return [Sequence.from_seqrecord(record) for record in (records if max_seqs == 0 else records[:max_seqs])]
        except Exception as e:
            raise ValueError(f"Error parsing sequence file: {e}")

    @staticmethod
    def _guess_format(filename: str) -> str:
        """Guess the sequence format from file extension"""
        ext = filename.lower().split('.')[-1]
        format_map = {
            'fa': 'fasta',
            'fasta': 'fasta',
            'sto': 'stockholm',
            'gb': 'genbank',
            'gbk': 'genbank',
            'phylip': 'phylip',
            'phy': 'phylip',
            'clustal': 'clustal',
            'aln': 'clustal',
            'embl': 'embl'
        }
        return format_map.get(ext, 'fasta')
