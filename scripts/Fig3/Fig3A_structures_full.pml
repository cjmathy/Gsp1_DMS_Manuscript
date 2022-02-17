########## Set working directory
cd /Users/cjmathy/gdrive/gsp1_dms/
delete all

########## Output image settings
set sphere_scale, 0.5
set cartoon_gap_cutoff, 0
set transparency, 0.5
set cartoon_transparency, 0.3
set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set dash_width, 4
set dash_gap, 0.5
set dash_color, green

###########################################################
### Figure 3A-l: GTP-bound structure highlighting Phe's ###
###########################################################

########## Load Gsp1-GTP structure
load "data/pdbs/pdbs_ran/3m1i.pdb"
remove solvent
util.cbaw 3m1i
util.cbay resn GTP
show sticks, 3m1i and not polymer
hide everything, 3m1i and resn Mg

########## Focus on active site, switches, and C-terminal linker
hide cartoon, 3m1i and resi 10-20+56-151+190-219
deselect

########## Color regions of the backbone
sele switch1, 3m1i and resi 30-47
color skyblue, switch1 and (name c,n,ca)

sele conformation_change, 3m1i and resi 172-208
color cyan, conformation_change and (name c,n,ca)

sele phenylalanines, 3m1i and resi 28+54+159+163
color 0xCD2027, phenylalanines

sele other_tox, 3m1i and resi 50+156
color salmon, other_tox and elem C

########## Show backbone, sidechain heavy atoms for residue network
select grp1_gtp, 3m1i and resi 28+30-37+39+48+50+54+152-154+156+159+163
show sticks, grp1_gtp
dist HB_gtp, (grp1_gtp),(grp1_gtp), mode=2
hide labels

########## Save image
set_view (\
    -0.598238885,    0.023611631,    0.800935507,\
    -0.008250237,   -0.999658585,    0.023293843,\
     0.801228166,    0.007331165,    0.598274052,\
     0.000000000,    0.000000000, -119.450042725,\
    61.026451111,   -2.196880341,   78.275840759,\
    72.883842468,  166.016189575,  -20.000000000 )
clip slab, 100
png figures/Fig3/Fig3A_GTP.png, dpi = 1200, 0, 0, -1

########## Remove the GTP-bound model
hide everything, 3m1i
delete HB_gtp

###########################################################
### Figure 3A-r: GDP-bound structure highlighting Phe's ###
###########################################################

########## Load Ran-GDP structure
load "data/pdbs/pdbs_ran/3gj0.pdb"
remove solvent
util.cbaw 3gj0
util.cbay resn GDP
show sticks, 3gj0 and not polymer
hide everything, 3gj0 and resn Mg

########## Focus on active site, switches, and C-terminal linker
hide cartoon, 3gj0 and resi 1-18+54-149+188-207
deselect

########## Color regions of the backbone
sele switch1, 3gj0 and resi 28-45
color skyblue, switch1 and (name c,n,ca)

sele conformation_change, 3gj0 and resi 170-206
color cyan, conformation_change and (name c,n,ca)
set cartoon_transparency, 0.3, 3gj0 and resi 177-189

sele phenylalanines, 3gj0 and resi 26+52+157+161
color 0xCD2027, phenylalanines

sele other_tox, 3gj0 and resi 48+154
color salmon, other_tox and elem C

########## Show backbone, sidechain heavy atoms for residue network
select grp1_gdp, 3gj0 and resi 26+28-35+37+46+48+52+150-152+154+157+161+181+183+186  # full
show sticks, grp1_gdp
hide cartoon, 3gj0 and resi 150-152
cartoon loop, 3gj0 and resi 150-152
show cartoon, 3gj0 and resi 150-152

dist HB_gdp, (grp1_gdp),(grp1_gdp), mode=2
hide labels

########## Save image
set_view (\
    -0.598238885,    0.023611631,    0.800935507,\
    -0.008250237,   -0.999658585,    0.023293843,\
     0.801228166,    0.007331165,    0.598274052,\
     0.000000000,    0.000000000, -119.450042725,\
    61.026451111,   -2.196880341,   78.275840759,\
    72.883842468,  166.016189575,  -20.000000000 )
clip slab, 100
png figures/Fig3/Fig3A_GDP.png, dpi = 1200, 0, 0, -1
