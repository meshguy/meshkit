###############################################################################
# resize.sh                                                                   #
#                                                                             #
# This program resizes a list of images for purpose of being viewed in        #
# browsers. The resized images replace the original images.                   #
#                                                                             #
# Author: Brett Rhodes                                                        #
# Date: 7/12/2013                                                             #
###############################################################################

#Max height of images allowed
MAX_H=480
#Max width of images allowed
MAX_W=960

#No file was given
if [ $# -lt 1 ]
then
  echo "Correct Usage is:"
  echo "  $0 file-names ..."
  exit 1
fi

#For each argument
for var in "$@"
do
  #Find it's height
  HEIGHT=`identify $var | sed -e 's/[^ ]* *[^ ]* *[0-9]*x\([0-9]*\).*/\1/'`
  #Find it's width
  WIDTH=`identify $var | sed -e 's/[^ ]* *[^ ]* *\([0-9]*\).*/\1/'`
  #Find the ratio to scale for Height
  HR=`expr $MAX_H \* 100 / $HEIGHT`
  #Find the ratio to scale for Width
  WR=`expr $MAX_W \* 100 / $WIDTH`

  #Determine which ratio is smaller
  if [ $HR -lt $WR ]
  then
    R=$HR
  else
    R=$WR
  fi

  #Resize, but do not enlarge image
  if [ $R -lt 100 ]
  then
    convert $var -resize $R% $var
  fi

done





