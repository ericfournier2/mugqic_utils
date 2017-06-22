#!/usr/bin/python

import sys,getopt
import csv
import os
import re

def main(argv):
   hub_name = 'hub_name'
   hub_short_label = 'hub_short_label'
   hub_long_label = 'hub_short_label'
   option_file = None
   input_dir = None
   try:
      opts, args = getopt.getopt(argv,"hn:s:l:o:",["name=", "shortlabel=", "longlabel=", "optionfile=", "help"])
   except getopt.GetoptError:
      print 'Usage: convert_to_hub.py [OPTION] <inputdir>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'Usage: convert_to_hub.py [OPTION] <inputdir>'
         print '-h, --help  '
         print '    Print this help.'
         print '-n VALUE, --name=VALUE'
         print '   Set the hub name.'
         print '-s VALUE, --shortlabel=VALUE'
         print '   Set the hub short label.'
         print '-l VALUE, --longlabel VALUE'
         print '   Set the hub log label.'
         print '-o VALUE, --optionfile=VALUE'
         print '   Give a file containing track option definitions. The file must have two'
         print '   columns, with the first column containing a regular expression, and the'
         print '   second column a track option to be set for tracks matching the regular'
         print '   expression.'
         sys.exit()
      elif opt in ("-n", "--name"):
         hub_name = arg
      elif opt in ("-s", "--shortlabel"):
         hub_short_label = arg
      elif opt in ("-l", "--longlabel"):
         hub_long_label = arg
      elif opt in ("-o", "--optionfile"):
         option_file = arg

   # Make sure we have an input directory.
   if len(args) == 0:
      print 'Usage: convert_to_hub.py [OPTION] <inputdir>'
      sys.exit(2)
   else:
       input_dir = args[0]


   # Generate the hub.txt file
   hub_content="""\
hub {hub_name}
shortLabel {hub_short_label}
longLabel {hub_long_label}
genomesFile genomes.txt
email Fournier.Eric.2@crchudequebec.ulaval.ca
descriptionUrl
""".format(**vars())

   with open(os.path.join(input_dir, "hub.txt"), "w") as w:
      w.write(hub_content)

   # Look for genome directories.
   genomes = next(os.walk(input_dir))[1]

   # Generate the genomes.txt file
   genome_content=""
   for genome in genomes:
      genome_content += "genome {0}\n".format(genome)
      genome_content += "trackDb {0}/trackDb.txt\n".format(genome)

   with open(os.path.join(input_dir, "genomes.txt"), "w") as w:
        w.write(genome_content)

   # Generate the trackDb files.
   # Read the option file.
   options = []
   if option_file != None:
      with open(option_file) as tsv:
         for line in csv.reader(tsv, dialect="excel-tab"):
            # Validate that the line has two elements.
            if len(line) == 2:
               option_pattern = re.compile(line[0])
               option_split = line[1].split(" ")
               if len(option_split) == 2:
                  options.append([option_pattern, option_split[0], option_split[1]])

   track_db = ""
   for genome in genomes:
      genome_dir = os.path.join(input_dir, genome)
      for track in next(os.walk(genome_dir))[2]:
         filename, fileext = os.path.splitext(os.path.basename(track))
         track_dict = {'shortLabel': filename, 'longLabel': filename, 'bigDataUrl': track}
         if (fileext == ".bw") or (fileext == ".bigwig"):
            track_dict['type'] = 'bigWig'
            track_dict['visibility'] = 'full'
         elif fileext == ".bb":
            track_dict['type'] = 'bigBed'
            track_dict['visibility'] = 'dense'
         elif fileext == ".bam":
            track_dict['type'] = 'bam'
            track_dict['visibility'] = 'dense'
         else:
            # Unrecognized file type. Skip it.
            continue

         for option in options:
            if option[0].search(track):
               track_dict[option[1]] = option[2]

         # Track must always appear first.
         track_db += "track {0}\n".format(filename)
         for key, value in track_dict.iteritems():
            track_db += "{0} {1}\n".format(key, value)
         track_db += "\n"

      with open(os.path.join(input_dir, genome, "trackDb.txt"), "w") as w:
         w.write(track_db)

if __name__ == "__main__":
    main(sys.argv[1:])