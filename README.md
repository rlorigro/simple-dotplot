# simple-dotplot
Very simple and configurable all-in-one dotplot program. No system dependencies other than make, CMake v3.10, and C++17. Create a dotplot directly from a pair of  FASTA/FASTQ files without having to run nucmer or any other alignment manually.

## Installation

```
git clone https://github.com/rlorigro/simple-dotplot.git
cd simple-dotplot
mkdir build
cd build
cmake ..
make -j [n_threads]
```

The executable is created in the `build` directory

## Usage 

```
Usage: ./dotplot [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -r,--ref TEXT REQUIRED      Path to FASTA/Q of ref sequences
  -q,--query TEXT REQUIRED    Path to FASTA/Q of query sequences
  -l,--min_length INT REQUIRED
                              Minimum length of match to be considered
  --mask_diagonal             Whether to zero all diagonal entries. 
                              Useful for self-dotplots, which have an extreme outlier on the diagonal
  -o,--output_dir TEXT REQUIRED
                              Path of directory to save output
```

Example:
```
./dotplot \
--ref /home/ryan/data/test/1.chm13.cenX.fasta \
--query /home/ryan/data/test/1.chm13.cenX.fasta \
--min_length 1024 \
--output_dir /home/ryan/data/test/test_simple_dotplot_mask_1024/ \
--mask_diagonal
```

## Example output

With `min_length=16`
![image](https://user-images.githubusercontent.com/28764332/124673796-5fcd7000-de6e-11eb-925d-c091b0e497ff.png)

With `min_length=1024`
![image](https://user-images.githubusercontent.com/28764332/124677724-e89bda00-de75-11eb-9e4e-3abb2b462ce9.png)


## Limitations

1. The directionality of your contigs must be set manually. The program will not try to search both the F and R complements.
2. At the moment there is no option to use MUMs instead of MEMs.

If these or any other limitations are causing you problems, please open an issue and I will update the repo when I get the chance.
