#!/usr/bin/env bash

# directory on the scratch file system where jobs will be run
dest="/home/anolan/scratch/thermal-structure/expr/00_CoupledInit"

# copy the source files 
rsync -upE initialize.* "${dest}/"

# Whole directories to copy to scratch
for dir in "run/" "sifs/" "params/"; do 
	rsync -rupE $dir $dest/$dir 
done

# folder structures to copy to scratch
for dir in "logs/" "result/"; do 
	rsync -ar --filter="-! */" $dir $dest/$dir
done 

# copy the mesh files 
rsync -upaR result/*/*/mesh.* $dest

for f in $(find ./ -name "*.coupled_init.time_profile"); do 
	# remove local path (./) from filename 
	file=${f/.\//} 
	
	# check if the symbolic link has already been created
	if [ ! -L "${dest}/${file}" ]; then
		# if not, create the symlink
		ln -s "${PWD}/${file}" "${dest}/${file}"
	fi

done