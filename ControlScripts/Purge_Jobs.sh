#! /bin/bash
# This script deletes the jobs from the gridsite
if [ "${1}"  == "--help" ] || [ "${1}" == "" ]; then
    echo "Options for running this script:"
    echo "./Purge_Jobs.sh --help                                        Prints this message"
    echo "./Purge_Jobs.sh --all                                         Purge ALL jobs"
    echo "./Purge_Jobs.sh --complete                                    Purge ONLY jobs which are completed"
    echo "./Purge_Jobs.sh --FORCEALL <host>[:tcpport]                   Removes all jobs at the specified site for the user"
    echo "./Purge_Jobs.sh --FORCEALL grid-ce.physik.rwth-aachen.de:8443 Removes all jobs at RWTH-Aachen for the user"
fi
if [ "${1}" == "--all" ]; then
    cat jobs_submitted | awk '{system("glite-ce-job-purge --noint " $1 " ;")}'
    cat jobs_complete | awk '{system("glite-ce-job-purge --noint " $1 " ;")}'
fi
if [ "${1}" == "--complete" ]; then
    cat jobs_complete | awk '{system("glite-ce-job-purge --noint " $1 " ;")}'
fi
if [ "${1}" == "--FORCEALL" ]; then
    glite-ce-job-purge --all --endpoint ${2}
fi
