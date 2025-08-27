from typing import Optional
from .color_scheme import ColorScheme

class SequenceFormatter:
    """Handles sequence formatting and coloring"""

    def __init__(self, color_scheme: ColorScheme):
        self.colors = color_scheme

    def _format_sequence_base(self, sequence: str, ncols: int, block_size: int, start_pos: int = 0, html: bool = False, consensus: Optional[str] = None, color_snps_only: bool = False) -> str:
        formatted = []
        ncols = 0 if start_pos >= len(sequence) else ncols
        sequence_slice = sequence[start_pos:] if ncols == 0 else sequence[start_pos:start_pos + ncols]
        consensus_slice = None
        if color_snps_only and consensus is not None:
            consensus_slice = consensus[start_pos:] if ncols == 0 else consensus[start_pos:start_pos + ncols]
            pairs = zip(sequence_slice, consensus_slice)
        else:
            pairs = zip(sequence_slice, [None]*len(sequence_slice))

        for i, (base, cons_base) in enumerate(pairs):
            if i > 0 and i % block_size == 0:
                formatted.append(' ')
            base_upper = base.upper()
            # Guard clause: skip coloring if SNP mode and matches consensus
            if color_snps_only and cons_base is not None and base_upper == cons_base.upper():
                formatted.append(base_upper)
                continue
            # Color if possible
            if html and base_upper in self.colors.html_colors:
                formatted.append(f'<span style="{self.colors.html_colors[base_upper]}">{base_upper}</span>')
                continue
            # Color for terminal output
            if not html and base_upper in self.colors.nucleotides:
                formatted.append(f"{self.colors.nucleotides[base_upper]}{base_upper}{self.colors.reset}")
                continue
            # No color, '-' or 'n'
            formatted.append(base_upper)
        return ''.join(formatted)

    def format_sequence(self, sequence: str, ncols: int, block_size: int, start_pos: int = 0, consensus: Optional[str] = None, color_snps_only: bool = False) -> str:
        return self._format_sequence_base(sequence, ncols, block_size, start_pos, html=False, consensus=consensus, color_snps_only=color_snps_only)

    def format_sequence_html(self, sequence: str, ncols: int, block_size: int, start_pos: int = 0, consensus: Optional[str] = None, color_snps_only: bool = False) -> str:
        return self._format_sequence_base(sequence, ncols, block_size, start_pos, html=True, consensus=consensus, color_snps_only=color_snps_only)

    def alignment_to_dataframe(self, sequences, sequence_names=None):
        """
        Converts a list of equal-length sequences to a pandas DataFrame.
        Each row is a sequence, each column is a position in the alignment.
        This is better for Plotly visualization with proper labels.

        Args:
            sequences: List of Sequence objects OR list of sequence strings
            sequence_names: Optional list of sequence names (used when sequences are strings)
        """
        import pandas as pd
        if not sequences:
            raise ValueError("No sequences provided.")

        # Handle both Sequence objects and strings
        if hasattr(sequences[0], 'sequence'):
            # Sequence objects
            length = len(sequences[0].sequence)
            for seq in sequences:
                if len(seq.sequence) != length:
                    raise ValueError("All sequences must be the same length.")

            # Create matrix
            matrix = [list(seq.sequence.upper()) for seq in sequences]

            # Create DataFrame with proper labels
            if sequence_names is None:
                sequence_names = [seq.header for seq in sequences]
        else:
            # String sequences
            length = len(sequences[0])
            for seq in sequences:
                if len(seq) != length:
                    raise ValueError("All sequences must be the same length.")

            # Create matrix
            matrix = [list(seq.upper()) for seq in sequences]

            # Create DataFrame with proper labels
            if sequence_names is None:
                sequence_names = [f"Seq_{i+1}" for i in range(len(sequences))]

        position_names = [i+1 for i in range(length)]

        df = pd.DataFrame(matrix, index=sequence_names, columns=position_names)
        return df
