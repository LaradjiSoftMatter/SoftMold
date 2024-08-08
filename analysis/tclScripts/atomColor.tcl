#index to index distances
#for example:
#indexDist 1 2 2 out.dat
proc indexDist { a b storeStep distFileName } {
	set n [molinfo top get numframes]
	set aAtom [atomselect top "index $a"]
	set bAtom [atomselect top "index $b"]
	set distFile [open $distFileName a]
	for { set i 0 } { $i < $n } { incr i } {
		$aAtom frame $i
		$bAtom frame $i
		set aCoord [lindex [$aAtom get {x y z}] 0 ]
		set bCoord [lindex [$bAtom get {x y z}] 0 ]
		set dist [distPBC $aCoord $bCoord [pbc get -now] ]
		set frameTime [expr $i*$storeStep]
		puts $distFile "$frameTime $dist"
	}
	close $distFile
}

#center of mass to center of mass distances
#for example:
#comDist [atomselect top "residue 1"] [atomselect top "residue 2"] 2 out.dat
proc comDist { aRes bRes storeStep distFileName } {
	set n [molinfo top get numframes]
	set distFile [open $distFileName a]
	for { set i 0 } { $i < $n } { incr i } {
		$aRes frame $i
		$bRes frame $i
		set aCom [center_of_mass $aRes]
		set bCom [center_of_mass $bRes]
		set dist [distPBC $aCom $bCom [pbc get -now] ]
		set frameTime [expr $i*$storeStep]
		puts $distFile "$frameTime $dist"
	}
	close $distFile
}


#radius of gyration of a molecule
#for example:
#gyrRadius [atomselect top "residue 1"] 2 out.dat
proc gyrRadius { sel storeStep gyrFileName } {
	set n [molinfo top get numframes]
	set gyrFile [open $gyrFileName a]
	# make sure this is a proper selection and has atoms
	if {[$sel num] <= 0} {
		error "gyrRadius: must have at least one atom in selection"
	}
	for { set i 0 } { $i < $n } { incr i } {
		$sel frame $i
		# gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
		set com [center_of_mass $sel]
		set sum 0
		foreach coord [$sel get {x y z}] {
			set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
		}
		set gyr [expr sqrt($sum / ([$sel num] + 0.0))]
		set frameTime [expr $i*$storeStep]
		puts $gyrFile "$frameTime $gyr"
	}
	close $gyrFile
}

#floating modulus operator
proc fmodulo {n m} {
	#puts "fmodulo: doing stuff"
	return [expr $n-$m*floor(1.0*$n/$m)]
}

#distance between points including periodic box coordinates (first 3 coordinates from pbcTools get)
proc distPBC {a b s} {
	set ax [lindex $a 0]
	set ay [lindex $a 1]
	set az [lindex $a 2]
	set bx [lindex $b 0]
	set by [lindex $b 1]
	set bz [lindex $b 2]
	set sx [lindex [lindex $s 0] 0]
	set sy [lindex [lindex $s 0] 1]
	set sz [lindex [lindex $s 0] 2]
	set fdx [expr $ax-$bx+($sx/2)]
	set fdy [expr $ay-$by+($sy/2)]
	set fdz [expr $az-$bz+($sz/2)]
	set gdx [fmodulo $fdx $sx]
	set gdy [fmodulo $fdy $sy]
	set gdz [fmodulo $fdz $sz]
	set dx [expr $gdx - ($sx/2)]
	set dy [expr $gdy - ($sy/2)]
	set dz [expr $gdz - ($sz/2)]
	return [expr sqrt($dx*$dx+$dy*$dy+$dz*$dz)]
}

#center of mass procedure from vmd tcl tutorial
proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}


proc massSet {mol index mass} {
	set sel [atomselect $mol "index $index"]
	$sel set mass
