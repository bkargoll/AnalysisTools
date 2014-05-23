#! /bin/bash
# This script deletes the jobs from the gridsite
if [ "${1}"  == "--help" ] || [ "${1}" == "" ]; then
    echo "Options for running this script:"
    echo "--help                    Prints this message"
    echo "--all                     Purge ALL jobs"
    echo "--complete                Purge ONLY jobs which are completed"
fi
if [ "${1}" == "--all" ]; then
    cat jobs_submitted | awk '{system("glite-ce-job-purge --noint " $1 " ;")}'
    cat jobs_complete | awk '{system("glite-ce-job-purge --noint " $1 " ;")}'
fi
if [ "${1}" == "--complete" ]; then
    cat jobs_complete | awk '{system("glite-ce-job-purge --noint " $1 " ;")}'
fi
