########## Set working directory
cd /Users/cjmathy/gdrive/gsp1_dms/
delete all

########## Output image settings
set sphere_scale, 0.5
set transparency, 0.25
set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set surface_quality, 3

########## Load Gsp1-GTP structure
load data/pdbs/pdbs_ran/3m1i.pdb, gsp1
remove solvent
color white, gsp1
util.cbaw (gsp1 and not polymer)
show sticks, gsp1 and not polymer
show spheres, gsp1 and resn Mg
hide everything, gsp1 and resi 183-220

########## Color different regions
color cyan, gsp1 and resi 172-182 # C-terminal linker
color slate, gsp1 and resi 96-171 # C-terminal lobe
sele hydrophobics, gsp1 and resi 88+89+90+92+94+100+103+106+110+117+119+120+121+138+140+146+148+165+166+167+170

########## Surface representation for hydrophobic side chains
create obj, hydrophobics
color slate, obj
show sticks, obj and not (name c,n)
show surface, obj
hide cartoon, obj

########## Plot
set_view (\
    -0.978801072,   -0.002088099,    0.204759374,\
    -0.201052070,   -0.179822475,   -0.962922990,\
     0.038829241,   -0.983685076,    0.175592378,\
    -0.000555622,    0.000037409, -121.200248718,\
    60.641456604,   -3.068675041,   73.981315613,\
    78.262344360,  164.117141724,  -20.000000000 )

png figures/FigS4/FigS4E_allosteric_lobe.png, dpi = 1200, 0, 0, -1
