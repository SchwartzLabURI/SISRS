#!/usr/bin/env python3

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

split_dir = post_processing_dir+"/Site_Splits"
good_split_dir = split_dir + "/Good_Splits"
good_split_lists = glob(good_split_dir+"/split*")

bad_split_dir = split_dir + "/Bad_Splits"
bad_split_lists = glob(bad_split_dir+"/bad_split*")

ref_annotation_dir = path.dirname(path.abspath(script_dir))+"/Reference_Genome/Annotations"
ref_annotation_files = glob(ref_annotation_dir+"/*.bed")

#Sort annotation BED files
for annotationFile in ref_annotation_files:
    sort_annotation = ['sort',
            '-k',
            '1,1',
            '-k2,2n',
            annotationFile,
            '-o',
            annotationFile]
    os.system(' '.join(sort_annotation))

sisrs_site_dir = post_processing_dir+"/SISRS_Sites"
sisrs_sites = sisrs_site_dir+"/"+ref_species+"_"+site_id+"_NoGaps_LocList.txt"

sort_sisrs = ['sort',
    '-u',
    '-o',
    '{}'.format(sisrs_sites),
    sisrs_sites]
subprocess.call(sort_sisrs)

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
                '-sorted',
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
                '|',
                'sort',
                '|'
                'uniq'
                '>',
                '{}'.format(output_anno)]
            os.system(' '.join(bed_command))
            print('Completed annotation transfer of '+annotation+' onto composite genome...')

            if(os.stat(output_anno).st_size == 0):
                count_file.write('Composite\tAll\tAll\t'+annotation+'\t0\n')
                outfile.write(annotation + " - No Sites\n")
            else:
                composite_annos.append(annotationFile)
                num_lines = sum(1 for line in open(output_anno))
                count_file.write('Composite\tAll\tAll\t'+annotation+'\t'+str(num_lines)+"\n")
                outfile.write(annotation + " - " + str(num_lines)+"\n")

sisrs_annos = []
with open(annotation_dir+"/Annotation_Counts.tsv","a+") as count_file:
    with open(post_log_dir+'/out_06_SISRS_Annotation',"w") as outfile:
        outfile.write("SISRS Annotation Counts:\n\n")
        for annotationFile in ref_annotation_files:
            annotation=path.basename(annotationFile).split('.')[0]
            input_anno ='{COMPOSITEANNODIR}/Composite_{ANNOTATION}_LocList.txt'.format(COMPOSITEANNODIR=composite_annotation_dir,ANNOTATION=annotation)
            output_anno ='{SISRSANNODIR}/SISRS_All_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,ANNOTATION=annotation)
            if(annotationFile in composite_annos):
                fetch_command = [
                    'comm',
                    '-12',
                    '{}'.format(input_anno),
                    '{}'.format(sisrs_sites),
                    '>',
                    '{}'.format(output_anno)]
                os.system(' '.join(fetch_command))
                if(os.stat(output_anno).st_size == 0):
                    count_file.write('SISRS\tAll\tAll\t'+annotation+'\t0\n')
                    outfile.write(annotation + " - No Sites\n")
                else:
                    sisrs_annos.append(annotationFile)
                    num_lines = sum(1 for line in open(output_anno))
                    count_file.write('SISRS\tAll\tAll\t'+annotation+'\t'+str(num_lines)+"\n")
                    outfile.write(annotation + " - " + str(num_lines)+"\n")
            else:
                count_file.write('SISRS\tAll\tAll\t'+annotation+'\t0\n')
                outfile.write(annotation + " - No Sites\n")

with open(annotation_dir+"/Annotation_Counts.tsv","a+") as count_file:
    with open(post_log_dir+'/out_07_SISRS_GoodSplit_Annotation',"w") as outfile:
        outfile.write("SISRS Good Split Annotation Counts:\n\n")
        for good_split in good_split_lists:
            split_num = path.basename(good_split).replace('split_','').replace('_LocList.txt','')
            for annotationFile in ref_annotation_files:
                annotation=path.basename(annotationFile).split('.')[0]
                input_anno ='{SISRSANNODIR}/SISRS_All_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,ANNOTATION=annotation)
                output_anno ='{SISRSANNODIR}/SISRS_Good_{SPLIT}_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,SPLIT=split_num,ANNOTATION=annotation)
                if(annotationFile in sisrs_annos):
                    fetch_command = [
                        'comm',
                        '-12',
                        '{}'.format(input_anno),
                        '{}'.format(good_split),
                        '>',
                        '{}'.format(output_anno)]
                    os.system(' '.join(fetch_command))
                    if(os.stat(output_anno).st_size == 0):
                        count_file.write('SISRS\tGood\t'+str(split_num)+'\t'+annotation+'\t0\n')
                        outfile.write('Good Split '+str(split_num)+" - " + annotation + " - No Sites\n")
                    else:
                        num_lines = sum(1 for line in open(output_anno))
                        count_file.write('SISRS\tGood\t'+str(split_num)+'\t'+annotation+'\t'+str(num_lines)+"\n")
                        outfile.write('Good Split '+str(split_num)+" - " + annotation + " - " + str(num_lines)+"\n")
                else:
                    count_file.write('SISRS\tGood\t'+str(split_num)+'\t'+annotation+'\t0\n')
                    outfile.write('Good Split '+str(split_num)+" - " + annotation + " - No Sites\n")

with open(annotation_dir+"/Annotation_Counts.tsv","a+") as count_file:
    with open(post_log_dir+'/out_08_SISRS_BadSplit_Annotation',"w") as outfile:
        outfile.write("SISRS Bad Split Annotation Counts:\n\n")
        for bad_split in bad_split_lists:
            split_num = path.basename(bad_split).replace('bad_split_','').replace('_LocList.txt','')
            for annotationFile in ref_annotation_files:
                annotation=path.basename(annotationFile).split('.')[0]
                input_anno ='{SISRSANNODIR}/SISRS_All_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,ANNOTATION=annotation)
                output_anno ='{SISRSANNODIR}/SISRS_Bad_{SPLIT}_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,SPLIT=split_num,ANNOTATION=annotation)
                if(annotationFile in sisrs_annos):
                    fetch_command = [
                        'comm',
                        '-12',
                        '{}'.format(input_anno),
                        '{}'.format(bad_split),
                        '>',
                        '{}'.format(output_anno)]
                    os.system(' '.join(fetch_command))
                    if(os.stat(output_anno).st_size == 0):
                        count_file.write('SISRS\tBad\t'+str(split_num)+'\t'+annotation+'\t0\n')
                        outfile.write('Bad Split '+str(split_num)+" - " + annotation + " - No Sites\n")
                    else:
                        num_lines = sum(1 for line in open(output_anno))
                        count_file.write('SISRS\tBad\t'+str(split_num)+'\t'+annotation+'\t'+str(num_lines)+"\n")
                        outfile.write('Bad Split '+str(split_num)+" - " + annotation + " - " + str(num_lines)+"\n")
                else:
                    count_file.write('SISRS\tBad\t'+str(split_num)+'\t'+annotation+'\t0\n')
                    outfile.write('Bad Split '+str(split_num)+" - " + annotation + " - No Sites\n")

with open(annotation_dir+"/Annotation_Counts.tsv","a+") as count_file:
    for bad_split in bad_split_lists:
        split_num = path.basename(bad_split).replace('bad_split_','').replace('_LocList.txt','')
        num_lines = sum(1 for line in open(bad_split))
        count_file.write('SISRS\tBad\t'+str(split_num)+'\tSISRS_Sites\t'+str(num_lines)+"\n")
    for good_split in good_split_lists:
        split_num = path.basename(good_split).replace('split_','').replace('_LocList.txt','')
        num_lines = sum(1 for line in open(good_split))
        count_file.write('SISRS\tGood\t'+str(split_num)+'\tSISRS_Sites\t'+str(num_lines)+"\n")

with open(annotation_dir+"/Annotation_Counts.tsv","a+") as count_file:
    for annotationFile in ref_annotation_files:
        annotation=path.basename(annotationFile).split('.')[0]
        ref_count_command = [
            'bedtools',
            'merge',
            '-i',
            '{}'.format(annotationFile),
            '|',
            'awk',
            "'{SUM",
            '+=',
            '$3-$2}',
            'END',
            '{print',
            "SUM}'",
            '-']
        ref_count = subprocess.check_output(' '.join(ref_count_command),shell=True).decode().strip()
        count_file.write('Reference\tAll\tAll\t'+annotation+'\t'+str(ref_count)+"\n")

os.system(' '.join(['cat']+good_split_lists+['| sort >',good_split_dir+'/Good_Split_LocList.txt']))
good_splits = good_split_dir+'/Good_Split_LocList.txt'

with open(annotation_dir+"/Annotation_Counts.tsv","a+") as count_file:
    for annotationFile in ref_annotation_files:
        annotation=path.basename(annotationFile).split('.')[0]
        input_anno ='{SISRSANNODIR}/SISRS_All_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,ANNOTATION=annotation)
        output_anno ='{SISRSANNODIR}/SISRS_AllGood_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,ANNOTATION=annotation)
        if(annotationFile in sisrs_annos):
            fetch_command = [
                'comm',
                '-12',
                '{}'.format(input_anno),
                '{}'.format(good_splits),
                '>',
                '{}'.format(output_anno)]
            os.system(' '.join(fetch_command))
            if(os.stat(output_anno).st_size == 0):
                count_file.write('SISRS\tGood\tAll\t'+annotation+'\t0\n')
            else:
                sisrs_annos.append(annotationFile)
                num_lines = sum(1 for line in open(output_anno))
                count_file.write('SISRS\tGood\tAll\t'+annotation+'\t'+str(num_lines)+"\n")
        else:
            count_file.write('SISRS\tGood\tAll\t'+annotation+'\t0\n')

bad_splits = bad_split_dir+'/Bad_Split_LocList.txt'

with open(annotation_dir+"/Annotation_Counts.tsv","a+") as count_file:
    for annotationFile in ref_annotation_files:
        annotation=path.basename(annotationFile).split('.')[0]
        input_anno ='{SISRSANNODIR}/SISRS_All_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,ANNOTATION=annotation)
        output_anno ='{SISRSANNODIR}/SISRS_AllBad_{ANNOTATION}_LocList.txt'.format(SISRSANNODIR=sisrs_annotation_dir,ANNOTATION=annotation)
        if(annotationFile in sisrs_annos):
            fetch_command = [
                'comm',
                '-12',
                '{}'.format(input_anno),
                '{}'.format(bad_splits),
                '>',
                '{}'.format(output_anno)]
            os.system(' '.join(fetch_command))
            if(os.stat(output_anno).st_size == 0):
                count_file.write('SISRS\tBad\tAll\t'+annotation+'\t0\n')
            else:
                sisrs_annos.append(annotationFile)
                num_lines = sum(1 for line in open(output_anno))
                count_file.write('SISRS\tBad\tAll\t'+annotation+'\t'+str(num_lines)+"\n")
        else:
            count_file.write('SISRS\tBad\tAll\t'+annotation+'\t0\n')
