set vmd_stride_convert_string {!#$%&'()*+,./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz|~}

#if {[info procs vmd_calculate_structure] == "vmd_calculate_structure"} {
#    rename vmd_calculate_structure vmd_calculate_structure_tmp
#    puts "did a rename."
#}
proc vmd_calculate_secstr {molid x i} {
    vmd_stride_full $molid $x $i
}
# If there was a rename, replace with the original
#if {[info procs vmd_calculate_structure_tmp] == \
#	"vmd_calculate_structure_tmp"} {
#    rename vmd_calculate_structure {}
#    rename vmd_calculate_structure_tmp vmd_calculate_structure
#}
#
proc luniq {L} {
    # removes duplicates without sorting the input list
    set t [list]
    foreach i $L {if {[lsearch -exact $t $i]==-1} {lappend t $i}}
    return $t
} ;# RS

set vmd_stride_firsttime 1
proc vmd_stride_full {molid x i} {
    global vmd_stride_firsttime
    if {$vmd_stride_firsttime} {
	set vmd_stride_firsttime 0
puts "In any publication of scientific results based in part or completely"
puts "on the use of the program STRIDE, please reference:"
puts "  Frishman,D & Argos,P. (1995) Knowledge-based secondary structure"
puts "  assignment.  Proteins:  structure, function and genetics, 23, 566-579."
    }

    global vmd_stride_convert_string env
    molinfo $molid get index
    set pfragsel [atomselect $molid "pfrag >= 0" ]
    $pfragsel set structure coil
    set pfrags [luniq [lsort -integer [$pfragsel get pfrag]]]
    set pnum [lindex $pfrags end]
    if {$pnum > 91} {
	puts "\nThere are $pnum protein fragments, but only 91 chain numbers"
	puts " are available.  Trying one protein fragment at a time"
	return [vmd_stride_simple $molid]
    }
    # are there fewer than 10,000 residues
    set resids [luniq [lsort [$pfragsel get residue]]]
    set num_resids [lindex $resids end]
    if {$num_resids > 9999} {
puts "\nThere are $num_resids protein in the molecule but STRIDE will only"
puts " analyze up to 10,000.  Trying one protein fragment at a time"
        return [vmd_stride_simple $molid]
    }

    # write the molecule to a PDB file
    global env
    set filename [pwd]/vmdstride0.in
    set dataname [pwd]/vmdstride0.out
    set outname  [pwd]/stride$x.dat

    vmd_stride_pdb_write $filename $pfrags $molid
    set arch_name [vmdinfo arch]
    set exec_name "$env(VMDDIR)/stride_$arch_name"
    catch {exec $exec_name -o $filename > exec0.out}
    exec awk {/^LOC/ {print $2,$4,$5,$7,$8}} exec0.out > $dataname
    #catch {exec $exec_name -o $filename | awk {/^LOC/ {print $2,$4,$5,$7,$8}} > $dataname 2> /dev/tty }
#    if [catch {exec $exec_name -o $filename | awk {/^LOC/ {print $2,$4,$5,$7,$8}} > $dataname 2> /dev/tty}] {
#	puts "First pass through STRIDE didn't work, trying one protein fragment at a time"
#	exec /bin/rm -f $filename $dataname
#	return [vmd_stride_simple $molid]
#    }
    set fileId [open $dataname r]
    set outId  [open $outname  a]

    while {[gets $fileId line] >= 0} {
	lassign $line ss resid1 chain1 resid2 chain2
	# convert from chain1 to pfrag
	set pfrag [string first $chain1 $vmd_stride_convert_string]
	set sel [atomselect $molid "pfrag $pfrag and resid $resid1 to $resid2"]
	set assign C
	switch -glob -- $ss {
	    AlphaHelix { set assign "H"}
	    310Helix { set assign "G"}
	    PiHelix { set assign "I"}
	    Strand { set assign "E"}
	    Bridge { set assign "B"}
	    Coil { set assign "C"}
	    Turn* { set assign "T"}
	    Gamma* { set assign "T"}
	    Disulfide {}
	    default { puts "unknown stride assignment $ss"}
	}
	$sel set structure [ldup [$sel num] $assign]
    puts $outId [format "%5d%1s%1s%1s%4d%1s%4d%1s%1s" $i " " $pfrag " " $resid1 " " $resid2 " " $assign]
    }
    close $fileId
    close $outId
    # get rid of the temp files
    #exec /bin/rm -f $filename $dataname
}

proc vmd_stride_pdb_write {filename pfrags molid} {
    global vmd_stride_convert_string
    # open the file (returns error if it couldn't)
    set fileId [open $filename w]
    puts $fileId [list REMARK STRIDE input file from VMD]

    # turn the string into an array
    set tmpstr $vmd_stride_convert_string
    for {set i 0} {$i < 91} {incr i} {
	set arr($i) [string range $vmd_stride_convert_string $i $i]
    }

    # write data from each pfrag
    foreach pfrag $pfrags {
	set sel [atomselect $molid "pfrag $pfrag"]
	set seldata [$sel get {
	 index name resname resid x y z occupancy beta segname}]
	set chain $arr($pfrag)
	foreach data $seldata {
	    lassign $data index name resname resid x y z \
		occupancy beta segname
	    incr index

   # the resname should end in column 3 (of the four) if < 4 characters
	    if {[string length $resname] < 4} {
		set resname "$resname "
	    }
	    puts $fileId [format "ATOM  %5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f %3s  %4s" $index $name " " $resname $chain $resid " " $x $y $z $occupancy $beta "   " $segname]
	}
    }
    puts $fileId END
    close $fileId
}


proc vmd_stride_simple {molid} {
    puts "ENTERING SIMPLE"
    # check that the molecule exists
    molinfo $molid get index

    # go along each pfrag
    set pfragsel [atomselect $molid "pfrag >= 0"]
    $pfragsel set structure coil
    set pfrags [luniq [lsort -integer [$pfragsel get pfrag]]]

    global env
    set filename [unique_file "$env(TMPDIR)/vmdstride0.in
    set dataname [unique_file "$env(TMPDIR)/vmdstride0.out
    foreach pfrag $pfrags {
	# get a pfrag
	set pfragsel [atomselect $molid "pfrag $pfrag"]
	# save to a file
	$pfragsel writepdb $filename

	set arch_name [vmdinfo arch]
	set exec_name  "$env(VMDDIR)/stride_$arch_name"
	if [catch {exec $exec_name -o $filename > $dataname 2> /dev/tty}] {
	    # get rid of the temp files
	    exec /bin/rm $filename $dataname
error "Could not execute STRIDE on data set! No secondary structure assignment"
	}

	# and read the LOC record inputs
	# 2th column for the residue, 4th for the start 7 for the end
	set fileId [open $dataname r]
	while {[gets $fileId line] >= 0} {
	    if [regexp "^LOC" $line] {
		lassign $line head ss name1 resid1 chain1 name2 resid2 chain2
		set sel [atomselect $molid "pfrag $pfrag and resid $resid1 to $resid2"]
                set assign C
		switch -glob -- $ss {
		    AlphaHelix { set assign "H"}
		    310Helix { set assign "G"}
		    PiHelix { set assign "I"}
		    Strand { set assign "E"}
		    Bridge { set assign "B"}
		    Coil { set assign "C"}
		    Turn* { set assign "T"}
		    Gamma* { set assign "T"}
		    default { puts "unknown stride assignment $ss"}
		}
#		puts [list $assign $resid1 $resid2]
		$sel set structure [ldup [$sel num] $assign]
	    }
	}
    }
    # and get rid of the temp files
    exec /bin/rm -f $filename $dataname
}


######## Not properly part of stride, but it is secondary structure related

## and a way to change the structure by hand
# eg: structure top "resid 5 to 12" helix
proc structure {molecule_id atom_selection structure_type} {
    set sel [atomselect $molecule_id $atom_selection]
    $sel set structure $structure_type
}

#mol delete [molinfo top]
for {set x 0} {$x <= 19} {incr x} {

    set a [expr 5*$x+1]
    set b [expr 5*$x+2]
    set c [expr 5*$x+3]
    set d [expr 5*$x+4]
    set e [expr 5*$x+5]

    mol new /site/raid3/student2/viki/dynamics/cab/step4_pbcsetup.psf
    mol addfile /site/raid3/student2/viki/dynamics/cab/traj/npt_$a.dcd waitfor all
    mol addfile /site/raid3/student2/viki/dynamics/cab/traj/npt_$b.dcd waitfor all
    mol addfile /site/raid3/student2/viki/dynamics/cab/traj/npt_$c.dcd waitfor all
    mol addfile /site/raid3/student2/viki/dynamics/cab/traj/npt_$d.dcd waitfor all
    mol addfile /site/raid3/student2/viki/dynamics/cab/traj/npt_$e.dcd waitfor all   
    set numframes [molinfo top get numframes]

    for {set i 0} {$i < $numframes} {incr i} {

        animate goto $i
        vmd_calculate_secstr [molinfo top] $x $i
    }
    mol delete [molinfo top]
}
