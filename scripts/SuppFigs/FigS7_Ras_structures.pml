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
set surface_quality, 1


load scripts/pymol_scripts/data2bfactor.py
load scripts/pymol_scripts/spectrumany.py
load scripts/pymol_scripts/spectrumbar.py

set_view (\
     1.000000000,    0.000000000,    0.000000000,\
     0.000000000,    1.000000000,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000000000,    0.000000000,  -34.027690887,\
     5.000000000,    0.000000000,    0.000000000,\
    26.827690125,   41.227691650,  -20.000000000 )

spectrumbar white,red, name=bar, ends=rounded
png figures/FigS7/FigS7_scale.png, dpi = 1200, 0, 0, -1

delete all

#############################################################################
# region of Ras allosteric switch (CaOAc2)
#    98-108 (Ras)
#    108-116 (Gsp1)
#
# residues in Ras allosteric network 1 from CaOAc2:
#    94+97+98+101+107+137 (Ras)
#    104+107+108+111+115+144 (Gsp1)
#
# residues in Ras allosteric network 2 from CaOAc2:
#    62+65+68+69+71+72+96+99+102 (Ras)
#    72+75+78+79+81+82+106+109+112 (Gsp1)
#
# residues with alternate dynamics in fraser et al 2011 pnas
#   20+61+62+68+72+86+88+92+94-96 (Ras)
#   29+71+72+78+82+96+98+102+104-106 (Gsp1)
#
# residues in fraser but not network 1 or 2
#   20+61+86+88+92+95 (Ras)
#   29+71+96+98+102+105 (Gsp1)
#
# Gorfe 2007 sites that show altered dynamics in presence of membrane
# (only GTPase fold sites, not including the C-terminal tail sites that
# are near the membrane)
#  103-109 (Ras site 1, allosteric switch)
#  113-117 (Gsp1 site 1, allosteric switch)   
#  117-126 (Ras site 2, in G4 region)
#  125-134 (Gsp1 site 2, in G4 region)
#  144-151 (Ras site 3, in G5 region)
#  151-158 (Gsp1 site 3, in G5 region)
#
#############################################################################

### NOTE for the Ras structures we are using the structures 2RGE from Buhrman et al
### (2007) Structure and 3K8Y from Buhrmen et al (2010) PNAS, as they are solved in the 
### same space group, and show the allosteric switch clearly. Other canonical Ras structures
### such as 5P21 and 1CTQ show the helix 3 to be closer to the "allosterically switched"
### structure, with the switch II loop more ordered, but notably these are in the space
### group P3221, in which the switch regions are modulated by crystal contacts (according
### to Buhrman et al (2007) Structure)

####### Load structures
load "data/pdbs/3k8y.pdb"
remove solvent
set_name 3k8y, HRas_allo
load "data/pdbs/2rge.pdb"
remove solvent
set_name 2rge, HRas_unsw
align HRas_allo, HRas_unsw
remove (HRas_unsw and not resi 58-76+92-111)

####### set main colors and style
util.cbaw HRas_allo
util.cbay HRas_allo and not polymer
show sticks, HRas_allo and not polymer
set_bond stick_radius, 0.5, HRas_allo and not polymer
show spheres, HRas_allo and resi 167-170
color palegreen, HRas_allo and resi 167+169
color yellow, HRas_allo and resi 168+170
util.cbaw HRas_allo and 
color palegreen, HRas_allo and resi 719 and symbol C

####### set unswitched region colors and style
util.cbaw HRas_unsw
set cartoon_transparency, 0.5, HRas_unsw
set stick_transparency, 0.5, HRas_unsw

####### select networks 1 and 2
sele buhrman_network1, resi 94+97+98+101+107+137
sele buhrman_network2, resi 62+65+68+69+71+72+96+99+102
sele fraser_addtl_res, resi 20+61+86+88+92+95
sele gorfe_site1, resi 103-109
sele gorfe_site2, resi 117-126
sele gorfe_site3, resi 144-151

sele show_sc, buhrman_network1 or buhrman_network2 or fraser_addtl_res
show sticks, (show_sc and (sidechain or name CA))
util.cbaw (show_sc and HRas_allo)
util.cbab (show_sc and HRas_unsw)

set_view (\
    -0.242301211,    0.353786528,   -0.903395593,\
     0.969384849,    0.050124198,   -0.240374655,\
    -0.039765224,   -0.933981776,   -0.355097443,\
    -0.000200240,    0.000167768, -106.284278870,\
    58.352512360,   -0.717430115,   74.399223328,\
    56.306423187,  156.306427002,  -20.000000000 )
clip slab, 100
png figures/FigS7/FigS7_switch.png, dpi = 1200, 0, 0, -1

remove HRas_unsw

sele show_sc, buhrman_network1 or buhrman_network2 or fraser_addtl_res or gorfe_site1 or gorfe_site2 or gorfe_site3
show sticks, (show_sc and (sidechain or name CA))
util.cbaw (show_sc and HRas_allo)


#### color sc shown positions by number of sensitive in Gsp1
alter show_sc, b=0
data2b_res HRas_allo, data/Ras_data/n_toxic_for_HRas_pos.txt
spectrum b, white red, show_sc, 0, 21
color red, show_sc and symbol O
color blue, show_sc and symbol N

########## Save images
set_view (\
    -0.324616909,    0.011721360,   -0.945772767,\
     0.873593271,   -0.379585385,   -0.304550588,\
    -0.362576157,   -0.925079644,    0.112983346,\
    -0.000180125,    0.000264883, -126.283668518,\
    61.124031067,   -2.442617416,   73.340156555,\
    76.283226013,  176.283279419,  -20.000000000 )
clip slab, 100
deselect
png figures/FigS7/FigS7_toxics.png, dpi = 1200, 0, 0, -1

#### color sc shown positions by number of sensitive in Gsp1
alter show_sc, b=0
data2b_res HRas_allo, data/Ras_data/n_activating_for_HRas_pos.txt
spectrum b, white red, show_sc, 0, 21
color red, show_sc and symbol O
color blue, show_sc and symbol N

png figures/FigS7/FigS7_activating.png, dpi = 1200, 0, 0, -1