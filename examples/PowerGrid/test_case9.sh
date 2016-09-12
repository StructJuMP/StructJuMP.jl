# runing test for case9 
echo "#####################################################"
echo "Teting (StructJuMP + Ipopt) Model: SC-ACOPF  Data: case9 with 2 contingencies"
julia scopf_structjump_main.jl Ipopt false  data/case9 2 8
#julia scopf_structjump_main.jl Ipopt true data/case9 2 8
echo "#####################################################"
echo "Teting (StructJuMP + PipsNlp) Model: SC-ACOPF Data: case9 with 2 contingencies"
julia scopf_structjump_main.jl PipsNlp false  data/case9 2 8
#julia scopf_structjump_main.jl PipsNlp true data/case9 2 8
echo "#####################################################"
echo "Teting (StructJuMP + PipsNlpSerial) Model: SC-ACOPF  Data: case9 with 2 contingencies"
julia scopf_structjump_main.jl PipsNlpSerial false  data/case9 2 8
#julia scopf_structjump_main.jl PipsNlpSerial true data/case9 2 8
echo "#####################################################"
echo "Teting (JuMP + Ipopt) Model: SC-ACOPF, Data: case9 with 2 contingencies"
julia scopf_main.jl data/case9 2 8
echo "#####################################################"
echo "Teting (JuMP + Ipopt) Model: ACOPF,  Data: case9"
julia acopf_main.jl data/case9 
echo "#####################################################"
