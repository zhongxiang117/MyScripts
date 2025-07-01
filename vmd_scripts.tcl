# VMD scripts
#
#> mkdir ~/.vmd ; cp PATH/TO/THIS_SCRIPT .
## add to ~/.vmdrc by `source THIS_SCRIPT`
#
## to use auto-completion, write `proc` names to
#> $env(VMDDIR)/vmd_completion.dat



proc xx_draw {args} {
    if {$args == "" || $args == "-h" || $args == "--help"} {
        puts "
Usage:

    link: https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node128.html

    * point {x y z}
    * line {x1 y1 z1} {x2 y2 z2} \[width 1.0] \[style <solid|dashed>]
    * cylinder {x1 y1 z1} {x2 y2 z2} \[radius 1.0] \[resolution 6] \[filled <yes|no>]
    * cone {basex basey basez} {tipz tipy tipz} \[radius 1.0] \[resolution 6]
    * triangle {x1 y1 z1} {x2 y2 z3} {x3 y3 z3}
    * trinorm {x1 y1 z1} {x2 y2 z3} {x3 y3 z3} {nx1 y1 z1} {nx2 ny2 nz3} {nx3 ny3 nz3}
    * tricolor {x1 y1 z1} {x2 y2 z3} {x3 y3 z3} {nx1 y1 z1} {nx2 ny2 nz3} {nx3 ny3 nz3} c1 c2 c3
    * sphere {x y z} \[radius 1.0] \[resolution 6]
    * text {x y z} ''text string'' \[size s] \[thickness t]
    * color colorId
    * color colorName
    * color trans_name
    * materials <on|off>
    * material name
    * delete id     : delete given id graphics primitive
    * delete all    : delete all graphics primitive
    * replace id    : next graphics primitive will replace given id
    * exists id     : check whether the primitive with the given id exists
    * list          : return a list of valid graphics ids
    * info id       : return command to recreate the graphics primitive with the given id

"
        return
    }
    draw $args
}



proc xx_tcl_gramma {} {
    puts "
Tcl/Tk Gramma

    for {init} {test} {increment} {commands}
    if {test} {commands}
    
    set file \[open ''file.txt'' w]
    puts \$file ''tasks''
    close \$file

    # check all environmental variables
    foreach key \[array names env] {puts \"\$key :  \$env(\$key)\"}


"
}



proc xx_draw_circle {{center {}} {vecx {1 0 0}} {vecy {0 1 0}} {radius 1.0} {resolution 36} {thickness 3} } {
    if {$center == "" || $center == "-h" || $center == "--help"} {
        puts "
Usage:

    >>> xx_draw_circle {center {0 0 0}} {vecx {1 0 0}} {vecy {0 1 0}} {radius 1.0} {resolution 36} {thickness 3}

    Spatially draw a circle on given center, two vectors, radius, resolution and thickness

    Idea:
        * two vectors define the planar where the circle will be drawn on by using segment lines
          - if they are unit vectors, radius will determine how big the circle is.
          - otherwise, radius will be acted as a factor for length between them with center point,
            thus a more general ellipse will be drawn on instead.
        * resolution defines the smooth of the circle
        * thickness is the segment line width
"
        return
    }
    set pi 3.1415926535897931
    for {set i 0} {$i < $resolution} {incr i} {
        set q0 [expr 2 * $pi * ($i-0.1) / $resolution]
        set t0 [vecadd $center [vecscale $vecx [expr $radius * cos($q0)]]]
        set r0 [vecadd $t0 [vecscale $vecy [expr $radius * sin($q0)]]]

        set q1 [expr 2 * $pi * ($i+1.1) / $resolution]
        set t1 [vecadd $center [vecscale $vecx [expr $radius * cos($q1)]]]
        set r1 [vecadd  $t1 [vecscale $vecy [expr $radius * sin($q1)]]]

        draw line $r0 $r1 width $thickness
    }
}



proc xx_draw_sphere_mesh {{center {}} {radius 1.0} {mesh 1} {resolution 36} {thickness 3} } {
    if {$center == "" || $center == "-h" || $center == "--help"} {
        puts "
Usage:

    >>> xx_draw_sphere_mesh {{center {0 0 0}} {radius 1.0} {mesh 1} {resolution 36} {thickness 3} }

    Spatially draw a sphere on given center, radius, resolution and thickness

    Idea:
        * mesh determines how many meshes will be used
        * resolution defines the smooth of the sphere
        * thickness is the segment line width
"
        return
    }
    set vecx [list $radius 0 0]
    set vecy [list 0 $radius 0]
    set vecz [list 0 0 $radius]

    xx_draw_circle $center $vecx $vecy $radius $resolution $thickness

    set pi 3.1415926535897931
    set dt [expr $pi / ($mesh-1)]
    for {set i 2} {$i<=$mesh} {incr i} {
        set x [expr $radius * cos($dt * ($i-1))]
        set y [expr $radius * sin($dt * ($i-1))]
        set vecm [list $x $y 0]
        xx_draw_circle $center $vecz $vecm $radius $resolution $thickness
    }
}



proc xx_draw_box {{point {}} {width 1} {length {}} {height {}}} {

    if {$point == "" || $point == "-h" || $point == "--help"} {
        puts "
Usage:

    >>> xx_draw_box {{point {0 0 0}} {width 1} {length {}} {height {}}}

    Idea:
        * box is drwan in the first quadrant of space (x+, y+, z+)
        * width in x+
        * length in y+, if not set, use width
        * height in z+, if not set, use width
        * thus:
          - if only point or width is given, cubic box is drawn;
          - otherwise, cuboid will be drawn instead

    #         zp________yz
    #          /       /|
    #         /       / | height
    #      xz/_______/pp|
    #        |*      |  /yp
    #   point|       | / width
    #        |_______|/
    #        xp      xy
    #          length

"
        return
    }

    if {$length eq ""} {set length $width}
    if {$height eq ""} {set height $width}

    set xb [lindex $point 0]
    set yb [lindex $point 1]
    set zb [lindex $point 2]

    set px [expr $xb + $width]
    set py [expr $yb + $length]
    set pz [expr $zb + $height]

    set xp [list $px $yb $zb]
    set xy [list $px $py $zb]
    set yp [list $xb $py $zb]
    set xz [list $px $yb $pz]
    set pp [list $px $py $pz]
    set zp [list $xb $yb $pz]
    set yz [list $xb $py $pz]

    draw line $point $xp
    draw line $point $yp
    draw line $point $zp
    draw line $xp $xy
    draw line $xp $xz
    draw line $xy $pp
    draw line $xy $yp
    draw line $pp $xz 
    draw line $pp $yz
    draw line $zp $yz
    draw line $zp $xz
    draw line $yz $yp
}



proc xx_draw_delete_ids {{start {}} {end {}}} {
    if {$start == "" || $start == "-h" || $start == "--help"} {
        puts "
Usage:

    >>> xx_draw_delete_ids start {end {}}

    delete draw objects by using their ids from start (inclusive) to end (inclusive)
    if end is not set, remove all objects later then start
"
        return
    }

    set ids [draw list]
    if {$end == ""} {
        set end [lindex $ids [expr [llength $ids] - 1]]
    }
    for {} {$start <= $end} {incr start} {draw delete $start}
}





