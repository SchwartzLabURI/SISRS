# *SISRS* v2.0 SNP Identification from Short Read Sequences


*SISRS* — pronounced “scissors” — is a program for identifying phylogenetically informative sites from next-generation whole-genome sequencing of multiple species. It identifies homologous sites without the need to do de novo assembly, annotation, and alignment. It identifies conserved regions by doing joint de novo assembly on multiple species. Sequencing reads are then aligned back to the contigs to identify variable sites.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

## Website

To learn how to use our software, and what software dependencies it has, we ask you to look at [SISRS Docs](https://schwartzlaburi.github.io/SISRS/).

## Support and Communication

If you have any questions about the software, please feel free to reach out to us on our github issues page @ [SISRS Github](https://github.com/SchwartzLabURI/SISRS/issues).

For other forms of communication we invite you to go to our lab's personal website @ [Schwartz Lab](https://schwartzlaburi.github.io/index.html).

## Tutorial

For sample data `unzip SISRS_Small.zip`

Assuming that you are running SISRS on an HPC with Slurm: `cd scripts`. Run all scripts from this folder.

`vim slurm_submissions/1_submit.slurm`

Change the modules, email address, and folder information as needed.

Repeat for remained slurm scripts
