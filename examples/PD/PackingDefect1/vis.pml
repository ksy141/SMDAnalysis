
load ../md.gro, struct

hide

cmd.select('PLhead', PLhead)
cmd.select('PLtail', PLtail)
cmd.select('TGhead', TGhead)
cmd.select('TGtail', TGtail)
deselect

show spheres, PLhead
color cPLhead, PLhead

show spheres, PLtail
color cPLtail, PLtail

show spheres, TGhead
color cTGhead, TGhead

show spheres, TGtail
color cTGtail, TGtail

load PLacyl.gro, PLacyl
color red, PLacyl

load TGacyl.gro, TGacyl
color orange, TGacyl

load TGglyc.gro, TGglyc
color purple, TGglyc

load Deep.gro, Deep
color pink, Deep

center struct
reset
util.performance(100)
rebuild

# ray 3600, 3600
# png 300ns.png


