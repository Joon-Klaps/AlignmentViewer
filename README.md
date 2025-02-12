# Alignment Viewer

A quick and simple FASTA alignment viewer for Python notebooks with colored sequence visualization.

## Features

- Color-coded nucleotide display
- Support for both Jupyter notebooks and terminal output
- Configurable display options (number of sequences, columns, ruler, etc.)
- Scrollable alignment view in notebooks
- Block-based sequence formatting

## Installation

```bash
pip install git+https://github.com/ECV-Lab-KULeuven/AlignmentViewer.git@dev
```

## Usage

```python
from AlignmentViewer import AlignmentViewer

# Display alignment from file
AlignmentViewer.display_alignment("data/alignment.fasta", nseqs=5, ncols=50)
```
![Example image](docs/example.png)

## Configuration

You can customize the display using the following parameters:

- `nseqs`: Number of sequences to display (default: 10)
- `ncols`: Number of columns to display (default: 100)
- `show_ruler`: Show position ruler (default: True)
- `block_size`: Size of sequence blocks (default: 10)
- `start_pos`: Starting position for display (default: 0)
- `container_height`: Height of the scrollable container in notebooks (default: "300px")
- `as_html`: Return HTML string instead of displaying in notebook (default: check)


## License

MIT License
