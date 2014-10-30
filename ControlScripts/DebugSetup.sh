#! /bin/bash
# This script deletes the jobs from the gridsite
if [ "${1}"  == "--help" ] || [ "${1}" == "" ]; then
    echo "#1 Recompiling your code with the debugging flag (note for externals or CommonUtils, you have to recomiple manually)"
    echo "cd Code; gmake clean; gmake CXXFLAGS+=-g "
    echo "#2 Use gdb to run the executable"
    echo "~/subs \"${PWD}/Code/Analysis.exe\" \"gdb ${PWD}/Code/Analysis.exe\" Set_1/Set_1.sh"
    echo '~/subs "$PWD/Code/Analysis.exe gdb" "$PWD/Code/Analysis.exe" Set_1/Set_1-GRID.sh'
fi
