axes location off
color Display Background iceblue
for {set i 0} {$i < 20} {incr i} {
    set filename snapshot_id$i.xyz
    mol new $filename autobonds off filebonds off
    scale by 2
    rotate x by -40
    rotate y by -40
    mol modstyle 0 $i CPK 1.0 0.5 30.0 30.0
    set filename snapshot_id$i.eps
    render PostScript $filename
    mol delete $i
}
quit
