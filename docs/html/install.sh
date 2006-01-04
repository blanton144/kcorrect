#!/bin/csh -f 

rsync -Larvz --delete $KCORRECT_DIR/docs/html/ ~/wwwblanton/kcorrect
