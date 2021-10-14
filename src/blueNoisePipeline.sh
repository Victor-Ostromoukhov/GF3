#!/bin/bash

# saner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o noclobber -o nounset

# -allow a command to fail with !’s side effect on errexit
# -use return value from ${PIPESTATUS[0]}, because ! hosed $?
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I’m sorry, `getopt --test` failed in this environment.'
    exit 1
fi


#Parameters
OPTIONS=d:n:o:
LONGOPTS=outputfolder:,optimtile:,vanillatile:,dim:,spp:,sobol_table:,dorender,prefix:


# -regarding ! and PIPESTATUS see above
# -temporarily store output to be able to check for errors
# -activate quoting/enhanced mode (e.g. by writing out “--options”)
# -pass arguments only via   -- "$@"   to separate them correctly
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

buildDir=cmake-build-release
pathDir=gk2path6d/bin
doInit=true
doOptim=true
doRender=false
sobol_dir=../data/sobol_init_tab.dat
outputFolder=.
begin=0 end=0
d=2 n=1 otile="" vtile="" prefix=tmp
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
	-o|--outputfolder)
	    outputFolder="$2"
	    shift 2
	    ;;
        -d|--dim)
            d="$2"
            shift 2
            ;;
        -n|--spp)
            n="$2"
            shift 2
            ;;
        --optimtile)
	    doInit=false
	    doOptim=false
            otile="$2"
            shift 2
            ;;
        --vanillatile)
	    doInit=false
            vtile="$2"
            shift 2
            ;;
        --prefix)
            prefix="$2"
            shift 2
            ;;
        --sobol_table)
            sobol_dir="$2"
            shift 2
            ;;
        --dorender)
	doRender=true
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

#exit 1

if $doInit; then
    echo "Initializing tile"
    echo "$buildDir/otDistanceMatrix -o $vtile --start $begin --end $end -n 128 --dirs $sobol_dir --dims 1 2 > $dist"
    $buildDir/otDistanceMatrix -o $vtile --start $begin --end $end -n 128 --dirs $sobol_dir --dims 1 2 > $dist
fi

if $doOptim; then
    echo "Optimizing tile"
    echo "$buildDir/Optimization/optimizeScreenSpacenD -i $vtile --loadDescDistancesFromFile $dist --dir_vectors $sobol_dir --dim $d -M $iterOptim -o $otile"
    $buildDir/Optimization/optimizeScreenSpacenD -i $vtile --loadDescDistancesFromFile $dist --dir_vectors $sobol_dir --dim $d -M $iterOptim -o $otile
fi

if $doRender; then
    echo "Rendering with tile"
    echo "$pathDir/path6d --prefix $prefix -s $n --sampler sobolpp -d $d --sobol_table $sobol_dir --seeds $otile --remap_dims "7 5:0 6:1" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell_sphere2_scene.txt --help"
    $pathDir/path6d --prefix $prefix -s $n --sampler sobolpp -d $d --sobol_table $sobol_dir --seeds $otile --remap_dims "7 5:0 6:1" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell_sphere2_scene.txt --help
fi

exit 0
















