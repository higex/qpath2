#!/usr/bin/env bash
#
# NAME
#
#    qpWsiToTiff - Convert a WSI to BigTIFF.
#
# SYNOPSIS
#    qpWsiToTIFF [<options>] <wsi_file> <out_file> [<temp_path>]
#
# DESCRIPTION
#    Convert a whole slide image (WSI) into a BigTIFF file
#    by extracting a single level and encoding it losslessly.
#    The layer can be specified via command line option and
#    cropped before being saved.
#
#    -h, --help
#        Output this documentation.
#    -l <n>, --level=<n>
#        Specify the level to be extracted (defaults to 0).
#    -c <x>,<y>,<width>,<height>, --crop=<x>,<y>,<width>,<height>
#        (x,y)-coordinates for the top-left corner of the crop region and its
#        width and height, in desired layer space.
#
#    wsi_file
#        A Whole Slide Image file understood by openslide.
#    out_file
#        Name of the result BigTIFF file.
#    temp_path
#        A folder for intermediate files. Defaults to /tmp but a file system
#        with faster access may be provided.
#
# DEPENDENCIES
#    This script relies on VIPS (https://github.com/jcupitt/libvips.git)
#    image processing library and OpenSlide (http://openslide.org) WSI I/O
#    library. VIPS should have been compiled with OpenSlide support.
#
# EXAMPLES
#    qpWsiToTiff input_wsi.ndpi output.tiff
#        Save the highest resolution (layer 0) from input_wsi.ndpi
#        in output.tiff
#
#    qpWsiToTiff --layer=3 -x 1000 -y 5000 -w 3000 -g 4000 input.ndpi out.tiff
#        Save the region from (1000,5000) to (4000,9000) from the 3rd layer
#        (downscale facto 2^3=8) of input.ndpi in out.tiff
#
# COPYRIGHT
#    Copyright (c) 2017 Vlad Popovici
#############################################################################

set -o errexit -o pipefail

# Exit codes from /usr/include/sysexits.h, as recommended by
# http://www.faqs.org/docs/abs/HTML/exitcodes.html
ex_ok=0           # Successful termination
ex_usage=64       # Command line usage error
ex_dataerr=65     # Data format error
ex_noinput=66     # Cannot open input
ex_nouser=67      # Addressee unknown
ex_nohost=68      # Host name unknown
ex_unavailable=69 # Service unavailable
ex_software=70    # Internal software error
ex_oserr=71       # System error (e.g., can't fork)
ex_osfile=72      # Critical OS file missing
ex_cantcreat=73   # Can't create (user) output file
ex_ioerr=74       # Input/output error
ex_tempfail=75    # Temp failure; user is invited to retry
ex_protocol=76    # Remote error in protocol
ex_noperm=77      # Permission denied
ex_config=78      # Configuration error

# All other errors
ex_unknown=1

# Base name of the sourcing script
script="`basename -- "$0"`"

# Print messages to standard error, with color if on an interactive
# terminal, then exits the sourcing script with the given exit code
# (default 1).
error() {
    printf '%s' "`basename -- "$0"`: " >&2

    if [ -t 1 ]
    then
        tput setf 4 || tput setaf 1
    fi

    # If the last parameter is a number, it's not part of the messages
    exit_code=1
    while true
    do
        if [ $# -eq 0 ]
        then
            break
        fi
        if [ $# -eq 1 ]
        then
            case "$1" in
                ''|*[!0-9]*)
                    ;;
                *)
                    exit_code="$1"
            esac
        fi
        printf '%s\n' "$1" >&2
        shift
    done

    if [ -t 1 ]
    then
       tput sgr0 # Reset formatting
    fi

    exit "$exit_code"
}


# Function to print documentation like the one above.
usage() {
    while IFS= read -r line || [ -n "$line" ]
    do
        case "$line" in
            '#!'*) # Shebang line
                ;;
            ''|'##'*|[!#]*) # End of comments
                exit "${1:-0}"
                ;;
            *) # Comment line
                printf '%s\n' "$line" >&2 # Remove comment prefix
                ;;
        esac
    done < "$0"
}


# Print messages to standard error, with color if on an interactive terminal.
warning() {
    if [ -t 1 ]
    then
        tput setf 4 || tput setaf 1
    fi

    printf '%s\n' "$@" >&2

    if [ -t 1 ]
    then
       tput sgr0 # Reset formatting
    fi
}

main() {
	# Process parameters
	declare -a parameters crop
	declare -i level x_crop y_crop w_crop h_crop crop
	while [ -n "${1:-}" ]
	do
		case $1 in
			-h | --help)
				usage
				exit 0
				;;
			-l | --level)
				level=$2
				shift 2
				;;
			--layer=*)
				layer=${1#*=}
				shift
				;;
			-c | --crop)
				crop=( $( echo $2 | sed 's/,/ /g' ) )
				shift 2
				;;
			--crop=*)
				crop=( $( echo ${1#*=} | sed 's/,/ /g' ) )
				shift
				;;
			--) # end of all options
				shift
				break
				;;
			-*) # illegal option
				error "ERROR: Unknown option $1"
				;;
			*)
				parameters+=( "$1" )
				shift
				;;
		esac
	done

	if [[ "${#parameters[@]}" -lt 2 ]]
	then
		warning "missing input or ouput file names"
		usage $ex_usage
	fi
	in_file=${parameters[0]}
	out_file=${parameters[1]}

	if [[ "${#parameters[@]}" -ge 3 ]]
	then
		tmp_path=${parameters[2]}
	else
		tmp_path="/tmp"
	fi

	# All set, do the work

	if [[ "${#crop[@]}" -gt 0 ]]
	then
		# if we crop, we go through an intermediate file
		if [[ "${#crop[@]}" -ne 4 ]]
		then
			error "incorrect crop region specification"
		fi
		vips openslideload --level ${level} \
			${in_file} ${tmp_path}/vips_tmp.v --vips-progress

		vips crop ${tmp_path}/vips_tmp.v \
			${out_file}[tile,tile_width=1024,tile_height=1024,compression=lzw] \
 			${crop[0]} ${crop[1]} ${crop[2]} ${crop[3]} --vips-progress 
		rm -f ${tmp_path}/vips_tmp.v
	else
		# no crop, save the final file directly
		vips openslideload --level ${level} \
			${in_file} \
			${out_file}[bigtiff,tile,tile_width=1024,tile_height=1024,compression=lzw] \
			--vips-progress
	fi
}


main "$@"

