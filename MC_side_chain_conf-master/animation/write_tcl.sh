codedir=C:/Users/79199/PycharmProjects/MC_side_chain_conf/animation
sed "s/YYY/$1/" $codedir/script_ex.tcl > $codedir/script_temp.tcl
sed "s/XXX/$2/" $codedir/script_temp.tcl >> $codedir/script.tcl