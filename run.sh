# source ./Apple_silcon.sh && cd utilities && make && cd .. && make 
source ./Apple_silcon.sh && cd utilities && make && cd ..
sudo bash -c "ulimit -s unlimited && ./utilities/ic.x >ic.log 2>&1"
make && ./main.x > main.log 2>&1 && ./utilities/density_power.x && cd utilities && ./dsp.x && ./fof.x && cd ..