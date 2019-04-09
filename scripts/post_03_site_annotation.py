import os
from os import path
import sys
from glob import glob
import subprocess
import stat
from itertools import islice

script_dir = sys.path[0]
processors = str(int(sys.argv[1]))
ref_species = sys.argv[2]
site_id = 'SISRS_Biallelic_NoMissing'

#Set directories based off of script folder location
post_processing_dir = path.dirname(path.abspath(script_dir))+"/Post_SISRS_Processing"

if(not path.isdir(post_processing_dir+"/Annotations")):
    os.mkdir(post_processing_dir+"/Annotations")
annotation_dir = post_processing_dir+"/Annotations"

if(not path.isdir(annotation_dir+"/Composite")):
    os.mkdir(annotation_dir+"/Composite")
composite_annotation_dir = annotation_dir+"/Composite"

if(not path.isdir(annotation_dir+"/SISRS")):
    os.mkdir(annotation_dir+"/SISRS")
sisrs_annotation_dir = annotation_dir+"/SISRS"

site_output_dir = post_processing_dir+"/SISRS_Sites"
good_split_dir = site_output_dir + "/Good_Splits"
bad_split_dir = site_output_dir + "/Bad_Splits"

ref_topology_dir = path.dirname(path.abspath(script_dir))+"/Reference_Topology"
ref_annotation_dir = path.dirname(path.abspath(script_dir))+"/Reference_Genome/Annotations"
ref_annotation_files = glob(ref_annotation_dir+"/*.bed")

sisrs_site_dir = post_processing_dir+"/SISRS_Sites"
composite_site_dir = post_processing_dir+"/Whole_Genome_Mapping"

post_log_dir = post_processing_dir+"/logFiles"

composite_annos = []

with open(annotation_dir+"/Annotation_Counts.tsv","w") as count_file:
    with open(post_log_dir+'/out_05_Composite_Annotation',"w") as outfile:
        outfile.write("Composite Genome Annotation Counts:\n\n")
        for annotationFile in ref_annotation_files:
            annotation=path.basename(annotationFile).split('.')[0]
            output_anno ='{COMPOSITEANNODIR}/Composite_{ANNOTATION}_LocList.txt'.format(COMPOSITEANNODIR=composite_annotation_dir,ANNOTATION=annotation)
            bed_command=[
                'bedtools',
                'intersect',
                '-a',
                '{}'.format(annotationFile),
                '-b',
                '{}'.format(composite_site_dir+"/WholeGenome_"+ref_species+"_Mapped_NonDup.bed"),
                '-wb',
                '|',
                'awk',
                "'{print",
                '$NF',
                "}'",
                '>',
                '{}'.format(output_anno)]
            os.system(' '.join(bed_command))
            if(os.stat(output_anno).st_size == 0):
                count_file.write('Composite\tAll\t'+annotation+'\t0\n')
                outfile.write(annotation + " - No Sites\n")
            else:
                composite_annos.append(annotationFile)
                num_lines = sum(1 for line in open(output_anno))
                count_file.write('Composite\tAll\t'+annotation+'\t'+str(num_lines)+"\n")
                outfile.write(annotation + " - " + str(num_lines)+"\n")
