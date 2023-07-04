#!/bin/bash

# Enforce running from ExaGO root dir
if [[ ! -f $PWD/buildsystem/build.sh ]]; then
	echo 'Please run this script from the top-level ExaGO source directory.'
    	exit 1
fi

# Need to know which platform we are configuring for
# Default to no_cluster
if [[ -z ${MY_CLUSTER} ]]
then
	echo "MY_CLUSTER is unset"
	echo "Try one of the following platforms: "
	echo $(ls -d ./buildsystem/spack/*/ | tr '\n' '\0' | xargs -0 -n 1 basename )
	exit 1	
fi

echo "$MY_CLUSTER" | awk '{print tolower($0)}'
# Use ${var,,} to convert to lower case
# There must be an existing folder for the cluster
if [ ! -d "./buildsystem/spack/${MY_CLUSTER}" ]
then
	echo "${MY_CLUSTER} did not match any directories in /buildsystem/spack/"
	echo "Try one of the following platforms: "
	echo $(ls -d ./buildsystem/spack/*/ | tr '\n' '\0' | xargs -0 -n 1 basename )
	exit 1
fi

# Make sure you are happy with your concretize before running this script
./buildsystem/spack/$MY_CLUSTER/install.sh

