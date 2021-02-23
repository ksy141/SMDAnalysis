
load md.gro, struct

hide

run sel.py
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


load PLacyl2.gro, PLacyl2
color red, PLacyl2

load TGacyl2.gro, TGacyl2
color orange, TGacyl2

load TGglyc2.gro, TGglyc2
color purple, TGglyc2

center struct
reset
util.performance(100)
rebuild

# ray 3600, 3600
# png 300ns.png

