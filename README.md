# LPSTyper  
***E. coli* LPS Outer Core Typing Tool**  

The tool identifies *E. coli* R1, R2, R3, R4, and K-12 outer core types in FASTA nucleotide sequences.

## Installation  

### Prerequisites  
Ensure you have the following installed:  

- **Python ≥ 3.9**  
  - Detailed installation instructions can be found [here](https://www.python.org/downloads/).  
- **Biopython** (Required for sequence processing)  
  - Install using:  
    ```
    pip install biopython
    ```

### Downloading LPSTyper  
You can obtain `LPSTyper.py` by either:  

- Downloading it manually from the repository  
- Using `git clone`:  
  ```
  git clone https://github.com/ellinium/LPSTyper.git
Usage
Run LPSTyper from the directory where it is located:

```
python LPSTyper.py <DIR>
```
where DIR is the path to a directory containing FASTA or FNA files (subdirectories are supported).

If running from another directory, provide the full path:
```
python /path/to/LPSTyper.py <DIR>
```
## Output
The results will be displayed on the screen.
Additionally, a CSV file with the results, LPSTyper_results.csv, will be saved in the provided directory DIR.
