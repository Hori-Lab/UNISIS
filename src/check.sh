cd ../check
#../build/bin/unisis input_check_force.toml
../build/bin/unisis input_T2S2HP.toml 1> ../src/out 2>&1
#../build/bin/unisis input_debug_ele.toml

cd ../src

cat out
