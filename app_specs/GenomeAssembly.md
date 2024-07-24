
# Application specification: GenomeAssembly

This is the application specification for service with identifier GenomeAssembly.

The backend script implementing the application is [App-GenomeAssembly.pl](../service-scripts/App-GenomeAssembly.pl).

The raw JSON file for this specification is [GenomeAssembly.json](GenomeAssembly.json).

This service performs the following task:   Assemble reads into a set of contigs

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_ids | SRR ID | string  |  |  |
| reference_assembly | Contig file | WS: Contigs  |  |  |
| recipe | Assembly recipe | enum  |  | auto |
| pipeline | Assembly pipeline arguments | string  |  |  |
| min_contig_len | Minimal output contig length | int  |  | 300 |
| min_contig_cov | Minimal output contig coverage | float  |  | 5 |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| debug | Activate debugging messages | int  |  |  |

