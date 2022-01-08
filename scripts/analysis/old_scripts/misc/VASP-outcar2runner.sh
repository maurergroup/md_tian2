#!/bin/bash
show_help() {
cat <<EOF
USAGE: VASP-outcar2runner.sh [-wWlq] [-c "first comment"] [-c "second comment"] [path-to-OUTCAR] > input.data
-c: add comments to the runner data file
-h: help, print this message
-l: only do last frame (last step of an optimization, for example)
-q: "quiet", do not autoprint details about calculation in comments
-w: energies and forces without vdW (default)
-W: energies and forces with vdW
[path-to-OUTCAR]: defaults to OUTCAR
EOF
}
with_vdw=0
only_last_frame=0
quiet=0
while getopts ":wWc:lqh" opt; do
    case "$opt" in
        w) with_vdw=0; ;;
        W) with_vdw=1; ;;
        c) [[ "$comment" = "" ]] && comment="$OPTARG" || comment="${comment}
$OPTARG"; ;;
        l) only_last_frame=1; ;;
        q) quiet=1; ;;
        h) show_help; exit; ;;
    esac
done
shift $((OPTIND-1))

if [[ "$1" = "" ]]; then
    working_dir=${USER}@${HOSTNAME%%[0-9]*}:$(pwd)
else
    working_dir=${USER}@${HOSTNAME%%[0-9]*}:$(dirname $(readlink -f $1))
fi

awk -v with_vdw=$with_vdw -v comment="$comment" -v working_dir="$working_dir" -v quiet="$quiet" '
BEGIN {
    BTOA=0.529177
    ATOB=1.0/BTOA
    EVTOHA=0.0367493
    split(comment,acomment,"\n")
}
NR==1 { sub(/^ */,""); sub(/ *$/,""); vaspversion=$0}
/NIONS =/ {nions=$NF}
/^ *GGA / {gga=$3}
/^ *ISPIN / {ispin=$3}
/^ *ENCUT / {encut=$3}
/^ *LASPH / {lasph=$3}
/^ *PREC / {prec=toupper(substr($3,1,1))}
/^ *ISMEAR / {ismear=$3; sigma=$6} #ismear includes a trailing semicolon
/^ *IVDW / {ivdw=$3}
/^ *VDW_SR / {vdwsr=$3}
/^ *TITEL / {pseudo=pseudo " " $4}
/^ *KPOINTS:/ {$1=""; kpoints=$0} #kpoints is the COMMENT LINE (first line) in the KPOINTS file
/NKPTS/ {nkpts=$4}
/VRHFIN/ {
    numtypes++
    typeT[numtypes]=substr($2,2, length($2)-2)
}
/ions per type =/ {
    for (i=5;i<=NF;i++){ntype[i-4]=$i}
    ntypes=NF-4
    if (ntypes != numtypes) { print "ERROR! ntypes != numtypes..." > "/dev/stderr"; exit 1; }
}
/direct lattice vectors/ {
    getline; l11=$1*ATOB; l12=$2*ATOB; l13=$3*ATOB; 
    getline; l21=$1*ATOB; l22=$2*ATOB; l23=$3*ATOB; 
    getline; l31=$1*ATOB; l32=$2*ATOB; l33=$3*ATOB; 
}
/POSITION         *TOTAL-FORCE/ {
    getline
    c=0
    for (t=1;t<=ntypes;t++) {
        for (i=1;i<=ntype[t];i++) {
            c++
            getline; a[c]=typeT[t]; 
            x[c]=$1*ATOB;y[c]=$2*ATOB;z[c]=$3*ATOB; #these are only given with four decimals in the VASP output
            if (with_vdw==1) {
                fx[c]=$4*EVTOHA*BTOA;fy[c]=$5*EVTOHA*BTOA;fz[c]=$6*EVTOHA*BTOA; #with vdW correction
            }
            else {
                fx[c]=$7*EVTOHA*BTOA;fy[c]=$8*EVTOHA*BTOA;fz[c]=$9*EVTOHA*BTOA; #without vdW correction
            }
        }
    }
}
/Hirshfeld charges/ {
    getline
    getline
    qsum=0
    for (i=1; i<=nions; i++) {
       getline
       q[i]=$3
       qsum+=q[i]
    }
}
/free  energy/ {
    if (with_vdw==1) energy=$5
    else energy=$10
    frame++
    print "begin"
    if (quiet==0) {
        print "comment VASP " vaspversion
        print "comment GGA = " gga "; ENCUT = " encut "; IVDW = " ivdw "; VDW_SR = " vdwsr "; ISPIN = " ispin "; ISMEAR = " ismear " SIGMA = " sigma "; PREC = " prec "; LASPH = " lasph "; frame " frame
        if (with_vdw==1) print "comment E and F *with* vdW contribution"
        else print "comment E and F *without* vdW contribution"
        print "comment PSEUDO " pseudo "; KPOINTS " kpoints "; NKPTS = " nkpts
        print "comment WD " working_dir
    }
    for (i=1; i<=length(acomment); i++) {
        print "comment " acomment[i]
    }
    printf "lattice %10.5f %10.5f %10.5f\n", l11, l12, l13
    printf "lattice %10.5f %10.5f %10.5f\n", l21, l22, l23
    printf "lattice %10.5f %10.5f %10.5f\n", l31, l32, l33
    for (i=1; i<=nions; i++) {
        printf "%5s %11.5f %11.5f %11.5f %2s %10.6f %10.6f %14.8f %14.8f %14.8f \n","atom ",x[i],y[i],z[i],a[i],q[i],0.0,fx[i],fy[i],fz[i]
    }
    printf "energy %16.9f\n", energy*0.0367493 #ev to Hartree
    printf "charge %16.9f\n", qsum
    print "end"
}


' ${1-OUTCAR} | {
    if [[ "$only_last_frame" = 1 ]]; then
        tac | sed '/^begin/q' | tac
    else
        cat
    fi
}
