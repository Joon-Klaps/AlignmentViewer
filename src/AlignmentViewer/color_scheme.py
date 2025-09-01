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
