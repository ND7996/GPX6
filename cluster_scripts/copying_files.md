COPYING EN FILES
scp -r nsekhar@192.168.133.2:/home/nsekhar/stepwise/MUT/step1/MOUSE_OLD/level0/R181S/rep{00..29} . 
scp -r nsekhar@192.168.133.2:/home/nsekhar/stepwise/MUT/step1/MOUSE_OLD/level0/R181S/traj{00..29} .
for folder in D148E E143S F139L F48Y G102S H144Q H177Q I24L K87T N107S P142S R181S R99C S47A S4R T178A T52A T60A Y104F; do      mkdir -p /home/hp/results/MOUSE/level1/$folder/replica000;     scp nsekhar@192.168.133.2:/home/nsekhar/stepwise/MUT/step1/MOUSE/level1/$folder/replica000/fep_000_1.000.dcd         nsekhar@192.168.133.2:/home/nsekhar/stepwise/MUT/step1/MOUSE/level1/$folder/replica000/fep_050_0.000.dcd         nsekhar@192.168.133.2:/home/nsekhar/stepwise/MUT/step1/MOUSE/level1/$folder/replica000/fep_025_0.500.dcd         /home/hp/results/MOUSE/level1/$folder/replica000/; done

