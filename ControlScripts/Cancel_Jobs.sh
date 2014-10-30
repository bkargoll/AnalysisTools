#! /bin/bash
# This script cancels the jobs on the gridsite
cat jobs_submitted | awk '{system("glite-ce-job-cancel --noint " $1 " ;")}'
