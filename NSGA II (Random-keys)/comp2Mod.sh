echo "---> CPLEX_Studio126"

cplex_path="/opt/ibm/ILOG/CPLEX_Studio126/"

# CPLEX

export PATH=$PATH:$cplex_path$"cplex/bin/x86-64_linux"

export LIBRARY_PATH=$LIBRARY_PATH:$cplex_path$"cplex/lib/x86-64_linux/static_pic"

export C_INCLUDE_PATH=$C_INCLUDE_PATH:$cplex_path$"cplex/include"
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$cplex_path$"cplex/include"

# CONCERT (CPP)

export LIBRARY_PATH=$LIBRARY_PATH:$cplex_path$"concert/lib/x86-64_linux/static_pic"

export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$cplex_path$"concert/include"

# op-L

export PATH=$PATH:$cplex_path$"opl/bin/x86-64_linux"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$cplex_path$"op-l/bin/x86-64_linux"

echo " --> Compiling"
rm *.out
rm *outputs.txt
rm *.ods
rm *.lp
rm salvaop-*

g++ -O2 -g -fmax-errors=10 -o OUTPUT.out NSGAII.cpp -DIL_STD -lilocplex -lconcert -lcplex -lm -lpthread -std=c++11 -fPIC

clear

echo " --> Running the model!"
echo " --> It's working!"

#G_SLICE=always-malloc G_DEBUG=gc-friendly  valgrind -v --track-origins=yes --tool=memcheck --leak-check=full --num-callers=40 --log-file=valgrind.log ./OUTPUT.out 1_6-1 0.4 40 8 0.1 0.3   .> salvaop-1_6_1.txt 
			 #inst cross pop it elite mutation	
			 


./OUTPUT.out 2_6-5_N_1_ 0.5 120 18 0 0.4  
./OUTPUT.out 2_6-4_N_1_ 0.5 120 18 0 0.4 
./OUTPUT.out 2_6-3_N_1_ 0.5 120 18 0 0.4 
./OUTPUT.out 2_6-2_N_1_ 0.5 120 18 0 0.4  
./OUTPUT.out 2_6-1_N_1_ 0.5 120 18 0 0.4 

./OUTPUT.out 2_8-5_N_1_ 0.5 120 18 0 0.4 
./OUTPUT.out 2_8-4_N_1_ 0.5 120 18 0 0.4 
./OUTPUT.out 2_8-3_N_1_ 0.5 120 18 0 0.4 
./OUTPUT.out 2_8-2_N_1_ 0.5 120 18 0 0.4  
./OUTPUT.out 2_8-1_N_1_ 0.5 120 18 0 0.4 

./OUTPUT.out 2_10-5_N_1_ 0.5 120 18 0 0.4 
./OUTPUT.out 2_10-4_N_1_ 0.5 120 18 0 0.4 
./OUTPUT.out 2_10-3_N_1_ 0.5 120 18 0 0.4 
./OUTPUT.out 2_10-2_N_1_ 0.5 120 18 0 0.4  
./OUTPUT.out 2_10-1_N_1_ 0.5 120 18 0 0.4

#sh mod1.sh &
#sh mod2.sh 
#sh mod3.sh &
#sh mod4.sh &
#sh mod5.sh &

#wait





mv *_frontier.txt ./Reg-Frontiers