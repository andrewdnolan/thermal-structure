#!/usr/bin/env bash

# directory on the scratch file system where jobs will be run
dest="/home/anolan/scratch/thermal-structure/expr/02_surge2steady/"

# copy the source files
rsync -uE surge2steady.* "${dest}/"

# Whole directories to copy to scratch
for dir in "run/" "sifs/" "params/"; do
	rsync -rE --delete $dir $dest/$dir
done

# folder structures to copy to scratch
for dir in "logs/" "result/"; do
	rsync -ar --no-g --no-p --filter="-! */" $dir $dest/$dir
done

# copy the mesh files
rsync -uaR --no-g --no-p result/*/*/mesh.* $dest
# copy whatever restrat files are included in the git repo
rsync -uaR --no-g --no-p result/*/*/*.result $dest


# make symbolic links to the time profile files
for f in $(find ./ -name "*.surge2steady.time_profile"); do
	# remove local path (./) from filename
	file=${f/.\//}

	# check if the symbolic link has already been created
	if [ ! -L "${dest}/${file}" ]; then
		# if not, create the symlink
		ln -s "${PWD}/${file}" "${dest}/${file}"
	fi

done
