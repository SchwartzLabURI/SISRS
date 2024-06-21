# *SISRS* v2.0 SNP Identification from Short Read Sequences


*SISRS* — pronounced “scissors” — is a program for identifying phylogenetically informative sites from next-generation whole-genome sequencing of multiple species. It identifies homologous sites without the need to do de novo assembly, annotation, and alignment. It identifies conserved regions by doing joint de novo assembly on multiple species. Sequencing reads are then aligned back to the contigs to identify variable sites.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation, either version 3 of the License, 
or (at your option) any later version.
You may use and modify this software so long as you acknowledge its authors,
list changes you have made, and continue to use the GPL3 license.

This program is distributed in the hope that it will be useful, but without any warranty; 
without even the implied warranty of merchantability or fitness for a particular purpose. 
See the GNU General Public License for more details.

## Tutorial

For sample data `unzip SISRS_Small.zip`

On an HPC using slurm, run the scripts in the slurm_example_scripts folder
in numerical order.
When using your own data ensure that the path to sisrs is correct in your scripts.

Note: to pick up an analysis from the middle but moving to a new folder 
(i.e. leaving the original analysis intact):
* run slurm script 01
* run 02 adding *--link* to the sisrs command to copy over prior results
* run 04 adding *--link [prior run folder]* 
* run 05 adding *--link [prior run folder]* 
* run 06 (array or not) adding *--link [prior run folder]* 
* run 10 adding *--link [prior run folder]* 


## Support and Communication

If you have any questions about the software, please feel free to reach out to us on our github issues page @ [SISRS Github](https://github.com/SchwartzLabURI/SISRS/issues).

For other forms of communication we invite you to go to our lab's personal website @ [Schwartz Lab](https://schwartzlaburi.github.io/).

