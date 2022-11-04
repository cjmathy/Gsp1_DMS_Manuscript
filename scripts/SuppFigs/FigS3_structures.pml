# set working directory
cd /Users/cjmathy/gdrive/gsp1_dms/

delete all

set sphere_scale, 0.5
set cartoon_gap_cutoff, 0

set transparency, 0.5
set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off

set surface_quality, 3
#set surface_quality, 1
load "data/pdbs/pdbs_ran/3m1i.pdb"
remove solvent
util.cbaw 3m1i
util.cbay resn GTP
show sticks, 3m1i and not polymer
set_bond stick_radius, 0.5, 3m1i and not polymer
show spheres, 3m1i and resn Mg
sele p_loop, 3m1i and resi 19-26
sele switch_1, 3m1i and resi 34-47
sele switch_2, 3m1i and resi 67-81
sele g_regions, 3m1i and resi 124-127+152-154
sele c_extension, 3m1i and resi 181-208
deselect
color firebrick, p_loop
color forest, switch_1
color yelloworange, switch_2
color skyblue, g_regions
color cyan, c_extension
create p_loop_obj, p_loop
show surface, p_loop_obj
show sticks, p_loop_obj and not (name c,n)
hide cartoon, p_loop_obj
create switch_1_obj, switch_1
show surface, switch_1_obj
show sticks, switch_1_obj and not (name c,n)
hide cartoon, switch_1_obj
create switch_2_obj, switch_2
show surface, switch_2_obj
show sticks, switch_2_obj and not (name c,n)
hide cartoon, switch_2_obj
create g_regions_obj, g_regions
show surface, g_regions_obj
show sticks, g_regions_obj and not (name c,n)
hide cartoon, g_regions_obj
create c_extension_obj, c_extension
show surface, c_extension_obj
show sticks, c_extension_obj and not (name c,n)
hide cartoon, c_extension_obj
set_view (\
    -0.211697459,    0.619699657,   -0.755737662,\
    -0.520752370,   -0.725874901,   -0.449333817,\
    -0.827033937,    0.298430383,    0.476379305,\
     0.000000000,    0.000000000, -170.828399658,\
    61.026451111,   -2.196880341,   78.275840759,\
   120.828392029,  220.828369141,  -20.000000000 )
clip slab, 100
png figures/FigS3/FigS3_GTP.png, dpi = 1200, 0, 0, -1

delete all

load "data/pdbs/pdbs_ran/3gj0.pdb"
remove solvent
util.cbaw 3gj0
util.cbay resn GDP
show sticks, 3gj0 and not polymer
set_bond stick_radius, 0.5, 3gj0 and not polymer
show spheres, 3gj0 and resn Mg
sele p_loop, 3gj0 and resi 17-24
sele switch_1, 3gj0 and resi 32-45
sele switch_2, 3gj0 and resi 65-79
sele g_regions, 3gj0 and resi 122-125+150-152
sele c_extension, 3gj0 and resi 179-206
deselect
color firebrick, p_loop
color forest, switch_1
color yelloworange, switch_2
color skyblue, g_regions
color cyan, c_extension
create p_loop_obj, p_loop
show surface, p_loop_obj
show sticks, p_loop_obj and not (name c,n)
hide cartoon, p_loop_obj
create switch_1_obj, switch_1
show surface, switch_1_obj
show sticks, switch_1_obj and not (name c,n)
hide cartoon, switch_1_obj
create switch_2_obj, switch_2
show surface, switch_2_obj
show sticks, switch_2_obj and not (name c,n)
hide cartoon, switch_2_obj
create g_regions_obj, g_regions
show surface, g_regions_obj
show sticks, g_regions_obj and not (name c,n)
hide cartoon, g_regions_obj
create c_extension_obj, c_extension
show surface, c_extension_obj
show sticks, c_extension_obj and not (name c,n)
hide cartoon, c_extension_obj
set_view (\
    -0.211697459,    0.619699657,   -0.755737662,\
    -0.520752370,   -0.725874901,   -0.449333817,\
    -0.827033937,    0.298430383,    0.476379305,\
     0.000000000,    0.000000000, -170.828399658,\
    61.026451111,   -2.196880341,   78.275840759,\
   120.828392029,  220.828369141,  -20.000000000 )
clip slab, 100
png figures/FigS3/FigS3_GDP.png, dpi = 1200, 0, 0, -1