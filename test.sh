if [ ! -e "mathdir" ]; then
	math -run 'Write["stdout", $InstallationDirectory]; Exit[]' -noprompt > mathdir
fi

mathdir=`cat mathdir`

mathdir=`echo $mathdir | sed "s/^\"\(.*\)\"$/\1/"`

echo $mathdir
