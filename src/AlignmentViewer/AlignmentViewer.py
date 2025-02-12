#!/usr/bin/env python3

from dataclasses import dataclass
from typing import Dict, List, Union, TextIO, Optional
from pathlib import Path
from IPython.display import HTML, display
from Bio import SeqIO

class ColorScheme:
    """Color configuration for sequence display"""
    nucleotides: Dict[str, str]
    reset: str
    html_colors: Dict[str, str]  # New HTML colors

    @classmethod
    def default(cls) -> 'ColorScheme':
        return cls(
            nucleotides={
                'A': '\033[48;2;76;114;165m\033[30m',   # blue
                'G': '\033[48;2;72;163;101m\033[30m',   # green
                'C': '\033[48;2;208;105;74m\033[30m',   # red
                'T': '\033[48;2;225;199;47m\033[30m',   # yellow
            },
            reset='\033[0m',
            html_colors={
                'A': 'background-color: #4c72a5; color: black;',
                'G': 'background-color: #48a365; color: black;',
                'C': 'background-color: #d0694a; color: black;',
                'T': 'background-color: #e1c72f; color: black;',
            }
        )

class DisplayConfig:
    """Configuration for alignment display.

    Attributes:
        nseqs (int): Number of sequences to display at once. Defaults to 10.
        ncols (int): Number of columns to display at once. Defaults to 100.
        show_ruler (bool): Whether to show position ruler above sequences. Defaults to True.
        block_size (int): Number of residues in each block before spacing. Defaults to 10.
        start_pos (int): Starting position for display (0-based). Defaults to 0.
        container_height (str): Height of the display container in CSS format. Defaults to "300px".
        as_html (Optional[bool]): Whether to return HTML format. If None, format is auto-detected.
            Defaults to None.
    """
    """Configuration for alignment display"""
    nseqs: int = 10
    ncols: int = 100
    show_ruler: bool = True
    block_size: int = 10
    start_pos: int = 0
    container_height: str = "300px"
    as_html: Optional[bool] = None

class Sequence:
    """Represents a single sequence in the alignment"""
    header: str
    sequence: str

    @classmethod
    def from_seqrecord(cls, record) -> 'Sequence':
        """Create a Sequence from a BioPython SeqRecord"""
        return cls(
            header=record.id,
            sequence=str(record.seq)
        )

class SequenceReader:
    """Handles sequence file parsing using BioPython"""

    @staticmethod
    def parse(source: Union[str, Path, TextIO], max_seqs: int) -> List[Sequence]:
        """Parse sequence format from various input sources using BioPython SeqIO"""
        try:
            # Try to determine format automatically for files
            if isinstance(source, (str, Path)):
                format = SequenceReader._guess_format(str(source))
                with open(source) as f:
                    records = list(SeqIO.parse(f, format))
            else:
                # For file handles, default to FASTA
                records = list(SeqIO.parse(source, "fasta"))

            return [Sequence.from_seqrecord(record) for record in records[:max_seqs]]
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
        return format_map.get(ext, 'fasta')  # Default to FASTA if unknown

class SequenceFormatter:
    """Handles sequence formatting and coloring"""

    def __init__(self, color_scheme: ColorScheme):
        self.colors = color_scheme

    def format_sequence(self, sequence: str, ncols: int, block_size: int, start_pos: int = 0) -> str:
        """Format and colorize a sequence"""
        formatted = []
        sequence_slice = sequence[start_pos:start_pos + ncols]

        for i, base in enumerate(sequence_slice):
            if i > 0 and i % block_size == 0:
                formatted.append(' ')
            base_upper = base.upper()
            if base_upper in self.colors.nucleotides:
                formatted.append(
                    f"{self.colors.nucleotides[base_upper]}"
                    f"{base_upper}{self.colors.reset}"
                )
            else:
                formatted.append(base_upper)
        return ''.join(formatted)

    def format_sequence_html(self, sequence: str, ncols: int, block_size: int, start_pos: int = 0) -> str:
        """Format sequence with HTML styling"""
        formatted = []
        sequence_slice = sequence[start_pos:start_pos + ncols]

        for i, base in enumerate(sequence_slice):
            if i > 0 and i % block_size == 0:
                formatted.append(' ')
            base_upper = base.upper()
            if base_upper in self.colors.html_colors:
                formatted.append(
                    f'<span style="{self.colors.html_colors[base_upper]}">'
                    f'{base_upper}</span>'
                )
            else:
                formatted.append(base_upper)
        return ''.join(formatted)

class AlignmentViewer:
    """Displays colored sequence alignments"""

    def __init__(self, color_scheme: Optional[ColorScheme] = None):
        self.color_scheme = color_scheme or ColorScheme.default()
        self.formatter = SequenceFormatter(self.color_scheme)

    @staticmethod
    def display_alignment(alignment: Union[str, Path, TextIO, List[Sequence]],
                         config: Optional[DisplayConfig] = None,
                         **kwargs) -> None:
        """Display a colored alignment with optional configuration"""
        viewer = AlignmentViewer()

        # Create base config
        base_config = config or DisplayConfig()

        # Update config with any provided kwargs
        for key, value in kwargs.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)
            else:
                raise ValueError(f"Unknown configuration parameter: {key}")

        # Parse sequences
        sequences = (alignment if isinstance(alignment, list)
                    else SequenceReader.parse(alignment, base_config.nseqs))

        # Calculate layout
        max_header_len = max(len(seq.header) for seq in sequences)

        # Use config's as_html if set, otherwise detect
        as_html = (base_config.as_html if base_config.as_html is not None
                      else viewer._check_notebook())

        if as_html:
            viewer._display_notebook(sequences, base_config, max_header_len)
        else:
            viewer._display_terminal(sequences, base_config, max_header_len)

    def _check_notebook(self) -> bool:
        """Check if we're running in a Jupyter notebook"""
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                return True
            return False
        except NameError:
            return False

    def _display_notebook(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment in a notebook with scrollable container"""
        html_parts = []

        # CSS for the container
        css = f"""
        <style>
            .alignment-container {{
                font-family: monospace;
                white-space: pre;
                overflow-x: auto;
                max-height: {config.container_height};
                background-color: white;
                padding: 10px;
                border: 1px solid #ddd;
            }}
            .sequence-row {{
                line-height: 1.5;
                margin: 2px 0;
            }}
        </style>
        """

        html_parts.append(css)
        html_parts.append('<div class="alignment-container">')

        # Add ruler if requested
        padding = ' ' * (max_header_len + 1)
        if config.show_ruler:
            ruler = self._create_ruler(padding, config.ncols, config.start_pos, config.block_size)
            html_parts.append(f"<div class='sequence-row'>{padding}{ruler}</div>")

        # Add sequences
        for sequence in sequences:
            colored_seq = self.formatter.format_sequence_html(
                sequence.sequence,
                config.ncols,
                config.block_size,
                config.start_pos
            )
            html_parts.append(
                f"<div class='sequence-row'>"
                f"{sequence.header:<{max_header_len}} {colored_seq}"
                f"</div>"
            )

        html_parts.append('</div>')

        # Display the HTML
        display(HTML(''.join(html_parts)))

    def _display_terminal(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment in terminal"""
        padding = ' ' * (max_header_len + 1)

        if config.show_ruler:
            print(f"{padding}{self._create_ruler(padding, config.ncols, config.start_pos, config.block_size)}")

        for sequence in sequences:
            colored_seq = self.formatter.format_sequence(
                sequence.sequence,
                config.ncols,
                config.block_size,
                config.start_pos
            )
            print(f"{sequence.header:<{max_header_len}} {colored_seq}")

    def _create_ruler(self, padding: str, width: int, start_pos: int = 0, step: int = 10) -> str:
        """Create a ruler string with column numbers and spaces"""
        ticks = ''
        numbers = []
        pos = start_pos

        # Create number line and tick marks simultaneously
        for i in range(width):
            if i % step == 0:
                num = str(pos + i)
                space = ' ' * (step - len(num))
                space += ' ' if i > 0 else ''
                numbers.append(num + space)
                ticks += '|'
                if i > 0:
                    ticks += ' '
            else:
                ticks += '-'

        return ''.join(numbers) + '\n' + padding + ticks