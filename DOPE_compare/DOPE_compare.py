from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller import soap_protein_od
import pylab
import os
import os.path
import sys

model1=sys.argv[1]
model2=sys.argv[2]
chain_id1=sys.argv[3]
chain_id2=sys.argv[4]
iD=sys.argv[5]


#~~~~~~~~~~~~~~~~~~~~~~ * Dual-template modelling, model evaluation and automatic loop refinement *~~~~~~~~~~~~~~~~~~~~~~~~#

# given two known homologous structures >30% iD, make:
# 1. Simple template model
# 2. Dual template model 
# 3. Loop refined model
# 4. Plot of best models from all 3 and their DOPE scores.
# 5. Renames best fitting models.

# Requires: pylab, os, sys, modeller
# usage: python2.7 DOPE_compare.py model1 model2 chainID1 chainID2 query_ID
# Comments with * added by me, all the rest from modeller example scripts 

#Notes on usage:
# - pdb_95.bin  pdb_95_new.pir  pdb_95.pir must be in dir above script #
# script should be run in working directory
# model1 and model2 are pdb file names but input without .pdb extension #
# model1 should be the model you want for the initial model #
# chain_id is chain e.g A,B,C in capitals #
# ID needs to be the same as id in the header and second line of .ali file e.g##
#>P1;test_protein
#sequence:test_protein:::::::0.00: 0.00
#MATLGTGMVSANDY*
#



def get_dopest_score(file_in, opt):

	"""* function to identify model with best DOPE score *"""

    with open(file_in, "r") as log_file:
        log_file=log_file.read()
        if opt == "M": # if not loop :
            x=log_file.split(">> Summary of successfully produced models:")[1]  
            x=x.split(">> Summary of successfully produced loop models:")[0]
            #print x
        elif opt == "L": # if loop
            x=log_file.split(">> Summary of successfully produced loop models:")[1] 
        elif opt == "ML":
            x=log_file.split(">> Summary of successfully produced models:")[1]
            x=x.split(">> Summary of successfully produced loop models:")[1]
        else:
             print "get_dopest_score: not an option - Use NL or L"
        
        
        z=filter(None, x.split("\n")[3:])

        D=[]
        P=[]

        for i in z:
            i=i.split()
            pdb_filename = i[0]
            dope_score = i[2]
            D.append(dope_score)
            P.append(pdb_filename)

        dopest_score=min(D, key=float) # find most negative DOPE #
        
  
        for line in zip(P,D):
            if str(dopest_score) == line[1]:
                best_pdb=line[0]              
        return best_pdb
        
        
        
#* ~~~~~~~~~~~~~~~~~~~~~~~~~~~Single template modelling~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *#


#* write out first log file *#

stout1 = sys.stdout 
sys.stdout= open(str(iD) + '_single.log', 'w')


# basic model #
env = environ()
log.verbose()  
aln = alignment(env)

mdl = model(env, file=str(model1), model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes=str(model1), atom_files=str(model1) + '.pdb')
aln.append(file=str(iD) + '.ali', align_codes=str(iD))
aln.align2d()
aln.write(file=str(iD) + 'single.ali', alignment_format='PIR')
aln.write(file=str(iD) + '-' + str(model1) + 'single.pap', alignment_format='PAP')

a1 = automodel(env, alnfile=str(iD) + 'single.ali',
              knowns=str(model1), sequence=str(iD),
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a1.starting_model = 1
a1.ending_model = 5  # n models to try
a1.make()


#* close first stdout link *#
sys.stdout.close()
sys.stdout=stout1
Best_single_model=get_dopest_score(str(iD) + '_single.log', "M")

#* Print best model info *# 
print "Best_single_model:" + str(Best_single_model) + " " + str(iD) + '_single_best_model.pdb'

#* Rename best model pdb file *#
os.rename(str(Best_single_model), str(iD) + '_single_best_model.pdb')  ## must run in same dir as output ##

#* evaluate best single model *#
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

#* read model file *#
mdl4 = complete_pdb(env, str(iD) + '_single_best_model.pdb')

#* Assess single model with DOPE *#
s = selection(mdl4)   # all atom selection
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=str(iD) + '_single_model.profile',
              normalize_profile=True, smoothing_window=15)



# directories for input atom files
env.io.atom_files_directory = './:../atom_files'

#* read single model file *#
mdl2 = complete_pdb(env, str(model1) + '.pdb', model_segment=('FIRST:A', 'LAST:A'))

#* evaluate single model *#
s = selection(mdl2)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=str(model1) + '.profile',
              normalize_profile=True, smoothing_window=15)





#* ~~~~~~~~~~~~~~~~~~~~~~~~~~~Multiple template modelling~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *#

#* write out second log file *#
stout2 = sys.stdout 
sys.stdout= open(str(iD) + '_multi.log', 'w')


#* Align 3D templates to each other *#
log.verbose()
env = environ()
env.io.atom_files_directory = './:../atom_files/'

aln = alignment(env)
for (code, chain) in (model1, chain_id1), (model2, chain_id2):
    mdl3 = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl3, atom_files=code, align_codes=code+chain)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                rr_file='$(LIB)/as1.sim.mat', overhang=30,
                gap_penalties_1d=(-450, -50),
                gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                dendrogram_file=str(iD) + '.tree',
                alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
                feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
                improve_alignment=True, fit=True, write_fit=write_fit,
                write_whole_pdb=whole, output='ALIGNMENT QUALITY')

aln.write(file=str(iD) + "_3dalign.pap", alignment_format='PAP')
aln.write(file=str(iD) + "_3dalign.ali", alignment_format='PIR')

aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
            rr_file='$(LIB)/as1.sim.mat', overhang=30,
            gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
            gap_gap_score=0, gap_residue_score=0, dendrogram_file=str(iD) + '.tree',
            alignment_type='progressive', feature_weights=[0]*6,
            improve_alignment=False, fit=False, write_fit=True,
            write_whole_pdb=False, output='QUALITY')
            
            
            
#* Align query sequence to the 3D template alignment *#
env.libs.topology.read(file='$(LIB)/top_heav.lib')

# Read aligned structure(s) #
aln = alignment(env)
aln.append(file=str(iD) + "_3dalign.ali", align_codes='all') ##** 3d structural alignment put here 
aln_block = len(aln)

# Read aligned sequence(s) #
aln.append(file=str(iD) + '.ali', align_codes=str(iD)) ###** original alignment input here

# Structure sensitive variable gap penalty sequence-sequence alignment:
aln.salign(output='', max_gap_length=20,
            gap_function=True,   # to use structure-dependent gap penalty
            alignment_type='PAIRWISE', align_block=aln_block,
            feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
            gap_penalties_1d=(-450, 0),
            gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
            similarity_flag=True)

aln.write(file=str(iD) + '-mult.ali', alignment_format='PIR')
aln.write(file=str(iD) + '-mult.pap', alignment_format='PAP')


# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']





#* ~~~~~~~~~~~~~~~~~~~~~~~~~~~Multiple template modelling with automatic loop refinement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *#



#* make loop models whilst making automatic loop refinement and get dope scores for both *#
a = dope_loopmodel(env,
              alnfile  = str(iD) + '-mult.ali',     # alignment filename
              knowns   = (str(model1) + str(chain_id1),str(model2)+str(chain_id2)), sequence = str(iD),assess_methods=assess.DOPE,
              loop_assess_methods=assess.DOPE)
              
a.starting_model= 1                 # index of the first model
a.ending_model  = 10  #5                 # index of the last model (determines how many models to calculate)
                                    
a.md_level = None                   # No refinement of model

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 8           # Last loop model
a.loop.md_level       = refine.slow # Loop model refinement level

a.make()                            # do comparative modeling

#* close second stdout *#
sys.stdout.close()
sys.stdout=stout2


#* Get best multi model and rename best model pdb file *#
Best_multi_model=get_dopest_score(str(iD) + '_multi.log', "M")
print "Best_multi_model:" + str(Best_multi_model) + " " + str(iD) + '_multi_best_model.pdb'
os.rename(str(Best_multi_model), str(iD) + '_multi_best_model.pdb')  # run in same dir as output #

#* get best loop model and rename best model pdb file*#
Best_loop_model=get_dopest_score(str(iD) + '_multi.log', "ML")
print "Best_loop_model:" + str(Best_loop_model) + " " + str(iD) + '_multi_best_loop_model.pdb'
os.rename(str(Best_loop_model), str(iD) + '_multi_best_loop_model.pdb')  # run in same dir as output #


#* evaluate multi model *#

# read model file
md0 = complete_pdb(env, str(iD) + '_multi_best_model.pdb')
# Assess all atoms with DOPE:
s = selection(md0)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=str(iD) + '_multi_model.profile',
              normalize_profile=True, smoothing_window=15)


#* evaluate loop multi model, get profile file *#

# read model file
md2 = complete_pdb(env, str(iD) + '_multi_best_loop_model.pdb')

#* Assess multi-loop-refined-model with DOPE *# 
s = selection(md2)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=str(iD) + '_multi_model_loop_refined.profile',
              normalize_profile=True, smoothing_window=15)


#* Plot profiles to compare single, multi and loop refined model fit based on DOPE scores *#

def r_enumerate(seq):
    """Enumerate a sequence in reverse order"""
    # Note that we don't use reversed() since Python 2.3 doesn't have it
    num = len(seq) - 1
    while num >= 0:
        yield num, seq[num]
        num -= 1

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    f = open(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for t, res in r_enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(t, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals


#* ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plot all DOPE scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *#

import modeller

k = modeller.environ()
y = modeller.alignment(k, file=str(iD) + '-mult.ali')

model = get_profile(str(iD) + '_multi_model.profile', y[str(iD)])  #   multi model 
template = get_profile(str(iD) + '_single_model.profile', y[str(iD)])    #  single model
loop_refined = get_profile(str(iD) + '_multi_model_loop_refined.profile', y[str(iD)])


#* Plot the template and model profiles in the same plot for comparison *#

pylab.figure(1, figsize=(10,6))
pylab.xlabel('Alignment position')
pylab.ylabel('DOPE per-residue score')

pylab.plot(model, color='red', linewidth=2, label='Multiple templates:' + str(model1) + " " + str(model2))
pylab.plot(template, color='green', linewidth=2, label='Basic model:' + str(model1))
pylab.plot(loop_refined, color='blue', linewidth=2, label='Loop refinement(DOPE)')

pylab.legend()
pylab.savefig(str(iD) + 'dope_profile.png', dpi=65)
