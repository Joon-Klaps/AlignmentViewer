from typing import Dict

class ColorScheme:
    """Color configuration for sequence display"""
    def __init__(self, nucleotides: Dict[str, str], reset: str, html_colors: Dict[str, str]):
        self.nucleotides = nucleotides
        self.reset = reset
        self.html_colors = html_colors

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

    def get_plotly_colors(self) -> Dict[str, str]:
        """Extract hex colors from color scheme for Plotly"""
        color_map = {}
        for nucleotide, html_style in self.html_colors.items():
            # Extract hex color from style string like "background-color: #4c72a5; color: black;"
            if 'background-color:' in html_style:
                hex_color = html_style.split('background-color:')[1].split(';')[0].strip()
                color_map[nucleotide] = hex_color

        # Add default colors for common nucleotides if not present
        defaults = {
            'A': '#4c72a5',  # blue
            'G': '#48a365',  # green
            'C': '#d0694a',  # red
            'T': '#e1c72f',  # yellow
            'U': '#e1c72f',  # yellow (same as T)
            '-': '#f0f0f0',  # light gray for gaps
            'N': '#cccccc'   # gray for unknown
        }

        for nuc, color in defaults.items():
            if nuc not in color_map:
                color_map[nuc] = color

        return color_map
