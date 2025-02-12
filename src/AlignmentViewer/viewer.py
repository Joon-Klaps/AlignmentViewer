#!/usr/bin/env python3

from dataclasses import dataclass
from typing import Dict, List, Union, TextIO, Optional
from pathlib import Path
from IPython.display import HTML, display

@dataclass
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

@dataclass
class DisplayConfig:
    """Configuration for alignment display"""
    nseqs: int = 10
    ncols: int = 100
    show_ruler: bool = True
    block_size: int = 10
    start_pos: int = 0  # New parameter for starting position
    container_height: str = "300px"  # New parameter for scroll container height

@dataclass
class Sequence:
    """Represents a single sequence in the alignment"""
    header: str
    sequence: str

class FastaReader:
    """Handles FASTA file parsing"""

    @staticmethod
    def parse(source: Union[str, Path, TextIO], max_seqs: int) -> List[Sequence]:
        """Parse FASTA format from various input sources"""
        if isinstance(source, (str, Path)):
            with open(source) as f:
                return FastaReader._parse_handle(f, max_seqs)
        return FastaReader._parse_handle(source, max_seqs)

    @staticmethod
    def _parse_handle(handle: TextIO, max_seqs: int) -> List[Sequence]:
        sequences = []
        current_header = ''
        current_seq = []

        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    sequences.append(Sequence(current_header, ''.join(current_seq)))
                    if len(sequences) >= max_seqs:
                        break
                current_header = line
                current_seq = []
            elif line:
                current_seq.append(line)

        if current_header and current_seq:
            sequences.append(Sequence(current_header, ''.join(current_seq)))

        return sequences[:max_seqs]

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
        self._is_notebook = self._check_notebook()

    def _check_notebook(self) -> bool:
        """Check if we're running in a Jupyter notebook"""
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                return True
            return False
        except NameError:
            return False

    @staticmethod
    def display_alignment(alignment: Union[str, Path, TextIO, List[Sequence]],
                         config: Optional[DisplayConfig] = None,
                         **kwargs) -> None:
        """Display a colored alignment with optional configuration"""
        # Create viewer instance for this call
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
                    else FastaReader.parse(alignment, base_config.nseqs))

        # Calculate layout
        max_header_len = max(len(seq.header) for seq in sequences)
        padding = ' ' * (max_header_len + 1)

        if viewer._is_notebook:
            viewer._display_notebook(sequences, base_config, max_header_len)
        else:
            viewer._display_terminal(sequences, base_config, max_header_len)

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