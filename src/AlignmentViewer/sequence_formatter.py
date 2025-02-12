from .color_scheme import ColorScheme

class SequenceFormatter:
    """Handles sequence formatting and coloring"""

    def __init__(self, color_scheme: ColorScheme):
        self.colors = color_scheme

    def _format_sequence_base(self, sequence: str, ncols: int, block_size: int, start_pos: int = 0, html: bool = False) -> str:
        formatted = []
        ncols = 0 if start_pos >= len(sequence) else ncols
        sequence_slice = sequence[start_pos:] if ncols == 0 else sequence[start_pos:start_pos + ncols]

        for i, base in enumerate(sequence_slice):
            if i > 0 and i % block_size == 0:
                formatted.append(' ')
            base_upper = base.upper()
            if base_upper in (self.colors.html_colors if html else self.colors.nucleotides):
                if html:
                    formatted.append(f'<span style="{self.colors.html_colors[base_upper]}">{base_upper}</span>')
                else:
                    formatted.append(f"{self.colors.nucleotides[base_upper]}{base_upper}{self.colors.reset}")
            else:
                formatted.append(base_upper)
        return ''.join(formatted)

    def format_sequence(self, sequence: str, ncols: int, block_size: int, start_pos: int = 0) -> str:
        return self._format_sequence_base(sequence, ncols, block_size, start_pos, html=False)

    def format_sequence_html(self, sequence: str, ncols: int, block_size: int, start_pos: int = 0) -> str:
        return self._format_sequence_base(sequence, ncols, block_size, start_pos, html=True)
