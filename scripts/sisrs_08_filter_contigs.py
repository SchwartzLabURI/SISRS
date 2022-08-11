import sys

#get file paths to alignment_pi_locs_m25.txt and contigs_SeqLength.tsv
#path_to_seq_lengths: Path to .tsv where contig lengths are located.
#                 Ex: "C://Users/caleb/OneDrive/Desktop/contigs_SeqLength.tsv"
#path_to_sites: Path to .txt file where names of individual sites on contigs are located
#                 Ex: "C://Users/caleb/OneDrive/Desktop/alignment_pi_locs_m25.txt"
#min_threshold: Removes contigs whos percentage of informative sites falls below a certain value.
#                 Ex: If you want to filter out contigs that are less than 25% informative, min_threshold
#                     should be 0.25
#path_to_output: Path to write filtered_data.tsv to. Must include the final backslash.
#                 Ex: "C://Users/caleb/OneDrive/Desktop/"


path_to_seq_lengths, path_to_sites, min_threshold = sys.argv[1], sys.argv[2], sys.argv[3]

#Example:
#path_to_seq_lengths, path_to_sites, min_threshold, path_to_output = "C://Users/caleb/OneDrive/Desktop/contigs_SeqLength2.tsv", "C://Users/caleb/OneDrive/Desktop/alignment_pi_locs_m25.txt", 0.1, "C://Users/caleb/OneDrive/Desktop/"

#path_to_sites = "C://Users/caleb/OneDrive/Desktop/alignment_pi_locs_m25.txt"
#path_to_seq_lengths = "C://Users/caleb/OneDrive/Desktop/contigs_SeqLength.tsv.txt.txt"

#read in informative site names
with open(path_to_sites) as file:
   sites = file.readlines()

#read in contig names and their overall lengths
with open(path_to_seq_lengths) as file:
    contigs = file.readlines()


#This function works similar to grep in R or bash, but it searches for an input string through
#a list rather than a file/directory. Strings in the list 'stringList' that contain the substring 'searchString'
#will be returned in list format.

grep = lambda searchString, stringList: [j for j in stringList if j.__contains__(searchString)]


contig_names = [i.split(",")[0] for i in contigs]
info_sites = [len((grep(i, sites))) for i in contig_names][2:]
possible_sites = [i.split(",")[1] for i in contigs]
possible_sites = [int(i) for i in possible_sites[1:]]

min_threshold = float(min_threshold)


for i in range(len(info_sites)):
    if(float(info_sites[i]) / float(possible_sites[i]) >= min_threshold):
        print(contig_names[i] + "   " + str(float(info_sites[i]) / float(possible_sites[i])))



