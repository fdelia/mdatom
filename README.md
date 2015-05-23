mdatom
======

Jede Änderung am Originalprogramm wurde mit einem "xxx" im Kommentar vermerkt.

Einstellung für das Coupling: ntt=
```
0: klassisch
1: Geschw. von T abhängig
2: nach 6.6.1
3: nach 6.6.2
```

params.inp
```
1: TITLE
2: NAT    NTXI    NTXO
3: BOX[1] BOX[2]  BOX[3]
4: NBOX[1]NBOX[2] NBOX[3] DX  IG  TEMPI
5: NTT    TEMPO   TAUT    BOLTZ
6: NSTLIM T       DT      
7: AMAS   EPSLJ   SIGLJ   RCUTF
8: NTPR   NTWX    NTWXM   NTPW
9: NGR    RCUTG
10:DTCOll GAMMA   TAUP    BETAT
```

Autoren: Timo Jung, Fabio D'Elia


## How to plot?
```
./mdatom exercises/6.6.1/params.inp exercises/6.6.1/coords.inp >> results.txt
python ../scripts/energy.py results.txt >> resultsFor.txt
``
then use R-Code
