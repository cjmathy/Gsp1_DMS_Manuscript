########## Set working directory
cd /Users/cjmathy/gdrive/gsp1_dms/
delete all

########## Output image settings
set sphere_scale, 0.5
set cartoon_gap_cutoff, 0
set transparency, 0.5
set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set surface_quality, 3
#set surface_quality, 1

########## Load Gsp1 structure
load "data/pdbs/pdbs_ran/3m1i.pdb"
remove solvent
hide everything, 3m1i


########################################################
### Figure 2B: Overlap between toxics and GTPase region UNION WITH CONTACTS NUCLEOTIDE
###
### toxic positions (from toxic_positions.csv)
###     19-26+28+32+35-37+43-46+50+54+56+59+67-72+75+77+78+80+81+93+96+97+99+101+125+127+132+133+152+153+155+156+157+159+160+163+178+182-184+186+187+191+192+202+206+211+220
### 
### residues in the GTPase region
###     19-26+34-47+67-81+124-127+152-154
### 
### residues that contact the nucleotide (from residues_contacting_nucleotide.txt))
###     20-27+37-41+43+44+67-71+124+125+127+128+152-154
########################################################

########## Copy a new Gsp1 model
create Gsp1, 3m1i and (resi 1-180) or resn GTP or resn Mg
show cartoon, Gsp1
util.cbaw Gsp1
util.cbay Gsp1 and not polymer
show sticks, Gsp1 and not polymer
show spheres, Gsp1 and resn Mg

########## Define toxics, the GTPase region residues, and the intersection of these two sets
sele toxics, Gsp1 and resi 19-26+28+32+35-37+43-46+50+54+56+59+67-72+75+77+78+80+81+93+96+97+99+101+125+127+132+133+152+153+155+156+157+159+160+163+178+182-184+186+187+191+192+202+206+211+220
sele g_region, Gsp1 and resi 19-26+34-47+67-81+124-127+152-154
sele contacting, Gsp1 and resi 20-27+37-41+43+44+67-71+124+125+127+128+152-154
sele union, toxics and (g_region or contacting)
deselect

########## Set colors
color 0xCD2027, toxics
color skyblue, union

########## Show sticks and surface for all toxic sidechains
create sc, toxics
show sticks, sc and not (name c,n)
show surface, sc
hide cartoon, sc

########## Save images
set_view (\
    -0.946020663,    0.081380866,   -0.313705742,\
     0.054834783,   -0.913809240,   -0.402420580,\
    -0.319416732,   -0.397904485,    0.860021412,\
     0.000499176,   -0.000059113, -127.357170105,\
    57.047744751,   -7.983699799,   77.965713501,\
    77.344955444,  177.344955444,  -20.000000000 )
clip slab, 100
png figures/Fig2/Fig2B_front.png, dpi = 1200, 0, 0, -1

set_view (\
    -0.433410257,   -0.285594881,    0.854740679,\
     0.008850582,   -0.949756980,   -0.312854588,\
     0.901148856,   -0.128031597,    0.414161593,\
    -0.000000000,    0.000000000, -148.618103027,\
    59.137012482,   -1.684617996,   74.924827576,\
   117.171646118,  180.064559937,  -20.000000000 )
png figures/Fig2/Fig2B_back.png, dpi = 1200, 0, 0, -1

########## Remove the Fig2A model to prepare for Fig2B
delete Gsp1 or sc or toxics or g_region or contacting or union





########################################################
### Figure 2C: Overlap between toxics and core positions
### 
### toxic positions defined in toxic_positions.csv
###     19-26+28+32+35-37+43-46+50+54+56+59+67-72+75+77+78+80+81+93+96+97+99+101+125+127+132+133+152+153+155+156+157+159+160+163+178+182-184+186+187+191+192+202+206+211+220
### 
### residues in core defined in burial.csv:
###     13-20+23-25+28+29+32+47+52+54+56+61+63+65+67+68+70+81+82+85+87-94+99+100+103+106+110+117-124+133+148+149+151-153+158+162-167+176
########################################################

########## Copy a new Gsp1 model
create Gsp1, 3m1i and (resi 1-180) or resn GTP or resn Mg
show cartoon, Gsp1
util.cbaw Gsp1
util.cbay Gsp1 and not polymer
show sticks, Gsp1 and not polymer
show spheres, Gsp1 and resn Mg

########## Define toxics, the core residues, and the intersection of these two sets
sele toxics, Gsp1 and resi 19-26+28+32+35-37+43-46+50+54+56+59+67-72+75+77+78+80+81+93+96+97+99+101+125+127+132+133+152+153+155+156+157+159+160+163+178+182-184+186+187+191+192+202+206+211+220
sele core, Gsp1 and resi 13-20+23-25+28+29+32+47+52+54+56+61+63+65+67+68+70+81+82+85+87-94+99+100+103+106+110+117-124+133+148+149+151-153+158+162-167+176
sele intersect, toxics and core
deselect

########## Set colors
color lightorange, core
color 0xCD2027, intersect

########## Show sticks and surface for all core sidechains
create sc, core
show sticks, sc and not (name c,n)
show surface, sc
hide cartoon, sc

########## Save images
set_view (\
    -0.946020663,    0.081380866,   -0.313705742,\
     0.054834783,   -0.913809240,   -0.402420580,\
    -0.319416732,   -0.397904485,    0.860021412,\
     0.000499176,   -0.000059113, -127.357170105,\
    57.047744751,   -7.983699799,   77.965713501,\
    77.344955444,  177.344955444,  -20.000000000 )
clip slab, 100
png figures/Fig2/Fig2C_front.png, dpi = 1200, 0, 0, -1

set_view (\
    -0.433410257,   -0.285594881,    0.854740679,\
     0.008850582,   -0.949756980,   -0.312854588,\
     0.901148856,   -0.128031597,    0.414161593,\
    -0.000000000,    0.000000000, -148.618103027,\
    59.137012482,   -1.684617996,   74.924827576,\
   117.171646118,  180.064559937,  -20.000000000 )
png figures/Fig2/Fig2C_back.png, dpi = 1200, 0, 0, -1

