########## Set working directory
cd /Users/cjmathy/gdrive/gsp1_dms/
delete all

########## Load script for setting variable based on outside dataset
# obtained from http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
load scripts/pymol_scripts/data2bfactor.py


########## Output image settings
bg white
set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set orthoscopic, 1   # so that spheres are scaled correctly


########## Load Gsp1-GTP structure
load "data/pdbs/pdbs_ran/3m1i.pdb"
remove solvent

color black, 3m1i
hide everything, 3m1i and polymer
show lines, 3m1i and polymer and (not resi 181-220) and ((name ca or name c or name n))
set_bond line_width, 5, 3m1i and polymer and ((name ca or name c or name n))
set sphere_scale, 0.25, resn Mg
util.cbay resn GTP
color yellow, resn Mg

########## Set b to the number of toxic mutations at each position
select residues, 3m1i and polymer
deselect
alter residues, b=0
data2b_res 3m1i, scripts/Fig4/toxics.txt
alter n. CA, vdw=b/10

########## Toxics in the canonical GTPase regions or contacting the nucleotide
show spheres, n. CA and resi 19-26+35-37+43-46+67-72+75+77+78+80+81+125+127+152+153
color white, n. CA and resi 19-26+35-37+43-46+67-72+75+77+78+80+81+125+127+152+153

########## Other toxics

show spheres, n. CA and resi 28+32+50+54+56+59+93+96+97+99+101+132+133+155+156+157+159+160+163+178 
color red, n. CA and resi 28+32+50+54+156+157+159+163
color purple, n. CA and resi 101+155
color teal, n. CA and resi 56+59+93+96+97+99+132+133+160+178

origin

create atom1, resi 96 and n. CA
translate [10,0,0], atom1
create atom2, atom1
create atom3, atom1

translate [0,4.5,0], atom2
translate [0,8,0], atom3

alter atom1 and n. CA, vdw=2
alter atom2 and n. CA, vdw=1.5
alter atom3 and n. CA, vdw=1


########## Plot

set_view (\
    -0.819974363,    0.553497910,   -0.145817712,\
    -0.566804826,   -0.820685625,    0.072103865,\
    -0.079770342,    0.141769245,    0.986672640,\
    -4.328055382,    2.729551077, -141.068206787,\
    61.239246368,   -2.122699738,   78.462867737,\
    94.463439941,  194.463409424,  -20.000000000 )

clip slab, 100
origin

set ray_trace_mode, 2

png figures/Fig4/Fig4_wire.png, dpi = 1200, 0, 0, -1

# In Adobe Illustrator, run [Object > Image Trace > Make and Expand], and then [Ungroup]