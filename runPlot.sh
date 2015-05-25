./mdatom exercises/6.6.1/params.inp > results/results1.txt
./mdatom exercises/6.6.2/params.inp exercises/6.6.2/coords.inp > results/results2.txt
./mdatom exercises/6.6.3/params.inp > results/results3.txt
python ../scripts/energy.py results/results1.txt > results/resultsFor1.txt
python ../scripts/energy.py results/results2.txt > results/resultsFor2.txt
python ../scripts/energy.py results/results3.txt > results/resultsFor3.txt