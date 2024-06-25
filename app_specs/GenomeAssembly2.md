
# Application specification: GenomeAssembly2

This is the application specification for service with identifier GenomeAssembly2.

The backend script implementing the application is [App-GenomeAssembly2.pl](../service-scripts/App-GenomeAssembly2.pl).

The raw JSON file for this specification is [GenomeAssembly2.json](GenomeAssembly2.json).

This service performs the following task:   Assemble reads into a set of contigs

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_ids | SRR ID | string  |  |  |
| recipe | Assembly recipe | enum  |  | auto |
| racon_iter | Racon iterations | int  |  | 2 |
| pilon_iter | Pilon iterations | int  |  | 2 |
| trim | trim_reads | boolean  |  | 0 |
| normalize | normalize_reads | boolean  |  | 0 |
| min_contig_len | Minimal output contig length | int  |  | 300 |
| min_contig_cov | Minimal output contig coverage | float  |  | 5 |
| genome_size | Genome Size | string  |  | 5M |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| debug | Debug level | int  |  | 0 |

