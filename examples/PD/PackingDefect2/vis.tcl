axes location off
display depthcue off
display shadows on
display ambientocclusion on
display projection orthographic

color Display Background white
display resize 1024 1024

### DEFECTS
mol representation lines

mol new PLacyl2.gro
mol delrep 0 top
mol addrep top

mol new TGglyc2.gro
mol delrep 0 top
mol color colorID 11
mol addrep top

mol new TGacyl2.gro
mol delrep 0 top
mol color colorID 3
mol addrep top

#mol new Deep.gro
#mol delrep 0 top
#mol color colorID 9
#mol addrep top

### STRUCTURE
mol new "../md.gro"

### PL
set    PLhead "N C11 C12 C13 C14 C15 C16 "
append PLhead "P O11 O12 O13 O14 "
append PLhead "C1 C11 C2 C21 C3 C31 "
append PLhead "O11 O12 O21 O22 O31 O32 "
append PLhead "HA HB HS HX HY "
append PLhead "HN1 HN2 HN3 "
append PLhead "H1 H2 H3 H4 H5 H6 "
append PLhead "HO2 HO3 HO4 HO5 HO6 "
append PLhead "O2 O3 O4 O5 O6 "


for {set i 1} {$i < 6} {incr i} {
append PLhead "H1${i}A H1${i}B H1${i}C "}

set    PLtail "H91 H101 "
for {set i 2} {$i < 19} {incr i} {
append PLtail "C2${i} H${i}R H${i}S H${i}T "
append PLtail "C3${i} H${i}X H${i}Y H${i}Z "}

puts "$PLhead"
puts "$PLtail"

set PLhead_sel [atomselect top "resname POPC DOPE SAPI and name $PLhead"]
set PLtail_sel [atomselect top "resname POPC DOPE SAPI and name $PLtail"]

$PLhead_sel set segname "PLhead"
$PLtail_sel set segname "PLtail"


### TG
set    TGhead "C1 C11 C2 C21 C3 C31 "
append TGhead "O11 O12 O21 O22 O31 O32 "
append TGhead "HA HB HS HX HY "

set    TGtail " "
for {set i 2} {$i < 19} {incr i} {
append TGtail "C1${i} H${i}A H${i}B H${i}C "
append TGtail "C2${i} H${i}R H${i}S H${i}T "
append TGtail "C3${i} H${i}X H${i}Y H${i}Z "}

set TGhead_sel [atomselect top "resname TRIO TRIV TRIN TRIP and name $TGhead"]
set TGtail_sel [atomselect top "resname TRIO TRIV TRIN TRIP and name $TGtail"]

$TGhead_sel set segname "TGhead"
$TGtail_sel set segname "TGtail"


### MATERIAL (make it brighter)
material add copy AOShiny
material change ambient Material23 0.150000 
mol material Material23
mol representation VdW 1.000000 42.000000
#mol representation QuickSurf 0.800000 0.500000 0.500000 3.000000


### COLOR
color change rgb 10 0.430000 0.710000 1.000000
color change rgb 0  0.000000 0.430000 0.860000
color change rgb 7  0.140000 1.000000 0.140000

color Segname PLhead 10
color Segname PLtail 0
color Segname TGhead 7
color Segname TGtail 17


### ADD ONE BY ONE
mol showrep top 0 0

mol color Segname
mol selection "resname POPC DOPE SAPI"
mol addrep top

mol color Segname
mol selection "resname TRIO TRIP TRIV TRIN and z > 140"
mol addrep top

scale by 3.5

