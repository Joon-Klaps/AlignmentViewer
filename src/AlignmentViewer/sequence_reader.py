from typing import List, Union, TextIO
from pathlib import Path
from Bio import SeqIO
from .sequence import Sequence

class SequenceReader:
    """Handles sequence file parsing using BioPython"""

    @staticmethod
    def parse(source: Union[str, Path, TextIO, List], max_seqs: int) -> List[Sequence]:
        """Parse sequence format from various input sources using BioPython SeqIO"""
        try:
            # Handle list input (could be SeqRecord objects or Sequence objects)
            if isinstance(source, list):
                return SequenceReader._parse_list(source, max_seqs)

            # Handle file paths and file objects
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
    def _parse_list(source_list: List, max_seqs: int) -> List[Sequence]:
        """Parse a list that could contain SeqRecord objects or Sequence objects"""
        if not source_list:
            return []

        # Check the first item to determine the type
        first_item = source_list[0]

        # Check if it's a BioPython SeqRecord (has 'seq' and 'id' attributes)
        if hasattr(first_item, 'seq') and hasattr(first_item, 'id'):
            # It's a list of SeqRecord objects
            max_seqs = 0 if max_seqs > len(source_list) else max_seqs
            records_to_process = source_list if max_seqs == 0 else source_list[:max_seqs]
            return [Sequence.from_seqrecord(record) for record in records_to_process]

        # Check if it's already a Sequence object
        elif hasattr(first_item, 'header') and hasattr(first_item, 'sequence'):
            # It's a list of Sequence objects
            max_seqs = 0 if max_seqs > len(source_list) else max_seqs
            return source_list if max_seqs == 0 else source_list[:max_seqs]

        else:
            raise ValueError(f"Unsupported list item type: {type(first_item)}. Expected BioPython SeqRecord or Sequence objects.")

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
