### Visualize the single snapshot
### Defects: Lines
### PL/TG: Material23 VdW

color Display Background white
display depthcue off
display shadows on
display ambientocclusion on
display resize 1024 1024
display nearclip set 0.010000
axes location Off
material add copy AOShiny
material change ambient Material23 0.150000

### PL HEAD/TAIL
set head "(resname POPC and name N C12 H12A H12B C13 H13A H13B H13C C14 H14A H14B H14C C15 H15A H15B H15C C11 H11A H11B \
P O13 O14 O12 O11 C1 HA HB C2 HS O21 C21 O22 C3 HX HY O31 C31 O32) or "

append head "(resname DOPE and name N HN1 HN2 HN3 C12 H12A H12B C11 H11A H11B \
P O13 O14 O11 O12 C1 HA HB C2 HS O21 C21 O22 C3 HX HY O31 C31 O32) or "

append head "(resname SAPI and name C12 H2 O2 HO2 C13 H3 O3 HO3 C14 H4 O4 HO4 C15 H5 O5 HO5 C16 H6 O6 HO6 C11 H1 \
P O13 O14 O12 O11 C1 HA HB C2 HS O21 C21 O22 C3 HX HY O31 C31 O32)"

set tail "(resname POPC DOPE SAPI) and name "
for {set i 2} {$i < 23} {incr i} {append tail "C2$i H${i}R H${i}S H${i}T C3$i H${i}X H${i}Y H${i}Z "}
append tail "H91 H101"

### TRIO HEAD/TAIL
set thead "resname TRIO and name C1 C2 C3 C11 C21 C31 O11 O21 O31 O12 O22 O32 HA HB HX HY HS"
set ttail "resname TRIO and name "
for {set i 2} {$i < 19} {incr i} { append ttail "C1${i} C2${i} C3${i} H${i}R H${i}S H${i}T H${i}A H${i}B H${i}C H${i}X H${i}Y H${i}Z "}
append ttail " "


### Simulation trajectory
mol new md.gro
mol showrep top 0 off

[atomselect top "$head"] set segname "HEAD"
[atomselect top "$tail"] set segname "TAIL"

color Segname HEAD 20
color Segname TAIL 21
color change rgb 20 0.430 0.710 1.00
color change rgb 21 0.000 0.430 0.86

[atomselect top "$thead"] set segname "THEAD"
[atomselect top "$ttail"] set segname "TTAIL"

color Segname THEAD 22
color Segname TTAIL 23
color change rgb 22 0.140 1.000 0.140
color change rgb 23 1.000 1.000 0.000

mol color Segname
mol representation VdW 1.000 42.000
mol material Material23
mol selection "resname POPC DOPE SAPI"
mol addrep top
mol selupdate 1 0 1

mol color Segname
mol representation VdW 1.000 42.000
mol material Material23
mol selection "resname TRIO"
mol addrep top
mol selupdate 2 0 1

### PL acyl defect
mol new PLacyl2.gro
mol showrep top 0 off
mol color colorID 1
mol representation Lines
mol selection "all"
mol addrep top

### Deep defect
mol new Deep2.gro
mol showrep top 0 off
mol color colorID 9
mol representation Lines
mol selection "all"
mol addrep top

### TG glycerol defect
mol new TGglyc2.gro
mol showrep top 0 off
mol color colorID 11
mol representation Lines
mol selection "all"
mol addrep top

### TG acyl defect
mol new TGacyl2.gro
mol showrep top 0 off
mol color colorID 3
mol representation Lines
mol selection "all"
mol addrep top

### Make the first trajectory to be top
mol top 0
display resetview

scale by 1.2
scale by 1.2
scale by 1.2

# render TachyonInternal 900ns.tga /usr/bin/open %s


