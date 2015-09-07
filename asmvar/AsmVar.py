"""
This is the main program of AsmVar. It's the most top interface of all the
AsmVar's tool sets.
"""
import sys
import time
import profile

#### My Own Module
import check
import commandline_options as cmdopts
import executor as exe

def checking():

    pass

    return

def genotype():

    gnt_opt = cmdopts.genotype()
	# Calculat the variants' genotype and output
    gnt_exe = exe.VariantsGenotype(gnt_opt.vcffile, 
                                   [f.split(':') for f in loadList(gnt_opt.bamfile)], 
                                   gnt_opt.ref_fasta_file,
                                   gnt_opt)

def loadList(file_name):

    I = open(file_name)
    file_list = [l.strip('\n') for l in I]
    I.close()

    return file_list


###############################################################################
if __name__ == '__main__':

    runner = {'genotype':genotype}

    if len(sys.argv) == 1 or (sys.argv[1] not in runner):
        print >> sys.stderr, '[Usage] python %s' % sys.argv[0]
        print >> sys.stderr, '\tOption:\n\t\tgenotype'
        sys.exit(1)

    command = sys.argv[1]
    runner[command]()
    #profile.run("genotype()") # Debug
    print >> sys.stderr, '** %s ALL DONE %s **\n' % (command, time.asctime())
    print >> sys.stderr, '**>> For the flowers bloom in the desert <<**'

