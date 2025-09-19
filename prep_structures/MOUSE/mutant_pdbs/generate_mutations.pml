# PyMOL script for GPX6 mutations
# Load this script in PyMOL: File -> Run Script -> select this file

load /home/hp/nayanika/github/GPX6/mouseWT/GPX6mousecys.pdb, original

# Level 1: T54Q
copy level1, original
alter resi 54 and level1, resn='GLN'
h_add level1
save mutant_pdbs/GPX6_level01.pdb, level1
delete level1

# Level 2: T54Q, C49U
copy level2, original
alter resi 54 and level2, resn='GLN'
alter resi 49 and level2, resn='SEC'
h_add level2
save mutant_pdbs/GPX6_level02.pdb, level2
delete level2

# Level 3: T54Q, C49U, I24L
copy level3, original
alter resi 54 and level3, resn='GLN'
alter resi 49 and level3, resn='SEC'
alter resi 24 and level3, resn='LEU'
h_add level3
save mutant_pdbs/GPX6_level03.pdb, level3
delete level3

# Level 4: T54Q, I24L, F137Y
copy level4, original
alter resi 54 and level4, resn='GLN'
alter resi 24 and level4, resn='LEU'
alter resi 137 and level4, resn='TYR'
h_add level4
save mutant_pdbs/GPX6_level04.pdb, level4
delete level4

# Level 5: T54Q, I24L, F137Y, S47A
copy level5, original
alter resi 54 and level5, resn='GLN'
alter resi 24 and level5, resn='LEU'
alter resi 137 and level5, resn='TYR'
alter resi 47 and level5, resn='ALA'
h_add level5
save mutant_pdbs/GPX6_level05.pdb, level5
delete level5

# Level 6: T54Q, I24L, F137Y, S47A, T50A
copy level6, original
alter resi 54 and level6, resn='GLN'
alter resi 24 and level6, resn='LEU'
alter resi 137 and level6, resn='TYR'
alter resi 47 and level6, resn='ALA'
alter resi 50 and level6, resn='ALA'
h_add level6
save mutant_pdbs/GPX6_level06.pdb, level6
delete level6

# Level 7: T54Q, I24L, F137Y, S47A, T50A, G74A
copy level7, original
alter resi 54 and level7, resn='GLN'
alter resi 24 and level7, resn='LEU'
alter resi 137 and level7, resn='TYR'
alter resi 47 and level7, resn='ALA'
alter resi 50 and level7, resn='ALA'
alter resi 74 and level7, resn='ALA'
h_add level7
save mutant_pdbs/GPX6_level07.pdb, level7
delete level7

# Level 8: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q
copy level8, original
alter resi 54 and level8, resn='GLN'
alter resi 24 and level8, resn='LEU'
alter resi 137 and level8, resn='TYR'
alter resi 47 and level8, resn='ALA'
alter resi 50 and level8, resn='ALA'
alter resi 74 and level8, resn='ALA'
alter resi 144 and level8, resn='GLN'
h_add level8
save mutant_pdbs/GPX6_level08.pdb, level8
delete level8

# Level 9: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C
copy level9, original
alter resi 54 and level9, resn='GLN'
alter resi 24 and level9, resn='LEU'
alter resi 137 and level9, resn='TYR'
alter resi 47 and level9, resn='ALA'
alter resi 50 and level9, resn='ALA'
alter resi 74 and level9, resn='ALA'
alter resi 144 and level9, resn='GLN'
alter resi 39 and level9, resn='CYS'
h_add level9
save mutant_pdbs/GPX6_level09.pdb, level9
delete level9

# Level 10: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q
copy level10, original
alter resi 54 and level10, resn='GLN'
alter resi 24 and level10, resn='LEU'
alter resi 137 and level10, resn='TYR'
alter resi 47 and level10, resn='ALA'
alter resi 50 and level10, resn='ALA'
alter resi 74 and level10, resn='ALA'
alter resi 144 and level10, resn='GLN'
alter resi 39 and level10, resn='CYS'
alter resi 177 and level10, resn='GLN'
h_add level10
save mutant_pdbs/GPX6_level10.pdb, level10
delete level10

# Level 11: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A
copy level11, original
alter resi 54 and level11, resn='GLN'
alter resi 24 and level11, resn='LEU'
alter resi 137 and level11, resn='TYR'
alter resi 47 and level11, resn='ALA'
alter resi 50 and level11, resn='ALA'
alter resi 74 and level11, resn='ALA'
alter resi 144 and level11, resn='GLN'
alter resi 39 and level11, resn='CYS'
alter resi 177 and level11, resn='GLN'
alter resi 179 and level11, resn='ALA'
h_add level11
save mutant_pdbs/GPX6_level11.pdb, level11
delete level11

# Level 12: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F
copy level12, original
alter resi 54 and level12, resn='GLN'
alter resi 24 and level12, resn='LEU'
alter resi 137 and level12, resn='TYR'
alter resi 47 and level12, resn='ALA'
alter resi 50 and level12, resn='ALA'
alter resi 74 and level12, resn='ALA'
alter resi 144 and level12, resn='GLN'
alter resi 39 and level12, resn='CYS'
alter resi 177 and level12, resn='GLN'
alter resi 179 and level12, resn='ALA'
alter resi 104 and level12, resn='PHE'
h_add level12
save mutant_pdbs/GPX6_level12.pdb, level12
delete level12

# Level 13: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F, S5R
copy level13, original
alter resi 54 and level13, resn='GLN'
alter resi 24 and level13, resn='LEU'
alter resi 137 and level13, resn='TYR'
alter resi 47 and level13, resn='ALA'
alter resi 50 and level13, resn='ALA'
alter resi 74 and level13, resn='ALA'
alter resi 144 and level13, resn='GLN'
alter resi 39 and level13, resn='CYS'
alter resi 177 and level13, resn='GLN'
alter resi 179 and level13, resn='ALA'
alter resi 104 and level13, resn='PHE'
alter resi 5 and level13, resn='ARG'
h_add level13
save mutant_pdbs/GPX6_level13.pdb, level13
delete level13

# Level 14: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F, S5R, T52A
copy level14, original
alter resi 54 and level14, resn='GLN'
alter resi 24 and level14, resn='LEU'
alter resi 137 and level14, resn='TYR'
alter resi 47 and level14, resn='ALA'
alter resi 50 and level14, resn='ALA'
alter resi 74 and level14, resn='ALA'
alter resi 144 and level14, resn='GLN'
alter resi 39 and level14, resn='CYS'
alter resi 177 and level14, resn='GLN'
alter resi 179 and level14, resn='ALA'
alter resi 104 and level14, resn='PHE'
alter resi 5 and level14, resn='ARG'
alter resi 52 and level14, resn='ALA'
h_add level14
save mutant_pdbs/GPX6_level14.pdb, level14
delete level14

# Level 15: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F, S5R, T52A, R87T
copy level15, original
alter resi 54 and level15, resn='GLN'
alter resi 24 and level15, resn='LEU'
alter resi 137 and level15, resn='TYR'
alter resi 47 and level15, resn='ALA'
alter resi 50 and level15, resn='ALA'
alter resi 74 and level15, resn='ALA'
alter resi 144 and level15, resn='GLN'
alter resi 39 and level15, resn='CYS'
alter resi 177 and level15, resn='GLN'
alter resi 179 and level15, resn='ALA'
alter resi 104 and level15, resn='PHE'
alter resi 5 and level15, resn='ARG'
alter resi 52 and level15, resn='ALA'
alter resi 87 and level15, resn='THR'
h_add level15
save mutant_pdbs/GPX6_level15.pdb, level15
delete level15

# Level 16: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F, S5R, T52A, R87T, Q102S
copy level16, original
alter resi 54 and level16, resn='GLN'
alter resi 24 and level16, resn='LEU'
alter resi 137 and level16, resn='TYR'
alter resi 47 and level16, resn='ALA'
alter resi 50 and level16, resn='ALA'
alter resi 74 and level16, resn='ALA'
alter resi 144 and level16, resn='GLN'
alter resi 39 and level16, resn='CYS'
alter resi 177 and level16, resn='GLN'
alter resi 179 and level16, resn='ALA'
alter resi 104 and level16, resn='PHE'
alter resi 5 and level16, resn='ARG'
alter resi 52 and level16, resn='ALA'
alter resi 87 and level16, resn='THR'
alter resi 102 and level16, resn='SER'
h_add level16
save mutant_pdbs/GPX6_level16.pdb, level16
delete level16

# Level 17: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F, S5R, T52A, R87T, Q102S, N107S
copy level17, original
alter resi 54 and level17, resn='GLN'
alter resi 24 and level17, resn='LEU'
alter resi 137 and level17, resn='TYR'
alter resi 47 and level17, resn='ALA'
alter resi 50 and level17, resn='ALA'
alter resi 74 and level17, resn='ALA'
alter resi 144 and level17, resn='GLN'
alter resi 39 and level17, resn='CYS'
alter resi 177 and level17, resn='GLN'
alter resi 179 and level17, resn='ALA'
alter resi 104 and level17, resn='PHE'
alter resi 5 and level17, resn='ARG'
alter resi 52 and level17, resn='ALA'
alter resi 87 and level17, resn='THR'
alter resi 102 and level17, resn='SER'
alter resi 107 and level17, resn='SER'
h_add level17
save mutant_pdbs/GPX6_level17.pdb, level17
delete level17

# Level 18: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F, S5R, T52A, R87T, Q102S, N107S, P142S
copy level18, original
alter resi 54 and level18, resn='GLN'
alter resi 24 and level18, resn='LEU'
alter resi 137 and level18, resn='TYR'
alter resi 47 and level18, resn='ALA'
alter resi 50 and level18, resn='ALA'
alter resi 74 and level18, resn='ALA'
alter resi 144 and level18, resn='GLN'
alter resi 39 and level18, resn='CYS'
alter resi 177 and level18, resn='GLN'
alter resi 179 and level18, resn='ALA'
alter resi 104 and level18, resn='PHE'
alter resi 5 and level18, resn='ARG'
alter resi 52 and level18, resn='ALA'
alter resi 87 and level18, resn='THR'
alter resi 102 and level18, resn='SER'
alter resi 107 and level18, resn='SER'
alter resi 142 and level18, resn='SER'
h_add level18
save mutant_pdbs/GPX6_level18.pdb, level18
delete level18

# Level 19: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F, S5R, T52A, R87T, Q102S, N107S, P142S, R181S
copy level19, original
alter resi 54 and level19, resn='GLN'
alter resi 24 and level19, resn='LEU'
alter resi 137 and level19, resn='TYR'
alter resi 47 and level19, resn='ALA'
alter resi 50 and level19, resn='ALA'
alter resi 74 and level19, resn='ALA'
alter resi 144 and level19, resn='GLN'
alter resi 39 and level19, resn='CYS'
alter resi 177 and level19, resn='GLN'
alter resi 179 and level19, resn='ALA'
alter resi 104 and level19, resn='PHE'
alter resi 5 and level19, resn='ARG'
alter resi 52 and level19, resn='ALA'
alter resi 87 and level19, resn='THR'
alter resi 102 and level19, resn='SER'
alter resi 107 and level19, resn='SER'
alter resi 142 and level19, resn='SER'
alter resi 181 and level19, resn='SER'
h_add level19
save mutant_pdbs/GPX6_level19.pdb, level19
delete level19

# Level 20: T54Q, I24L, F137Y, S47A, T50A, G74A, H144Q, R39C, H177Q, T179A, Y104F, S5R, T52A, R87T, Q102S, N107S, P142S, R181S, F80Y
copy level20, original
alter resi 54 and level20, resn='GLN'
alter resi 24 and level20, resn='LEU'
alter resi 137 and level20, resn='TYR'
alter resi 47 and level20, resn='ALA'
alter resi 50 and level20, resn='ALA'
alter resi 74 and level20, resn='ALA'
alter resi 144 and level20, resn='GLN'
alter resi 39 and level20, resn='CYS'
alter resi 177 and level20, resn='GLN'
alter resi 179 and level20, resn='ALA'
alter resi 104 and level20, resn='PHE'
alter resi 5 and level20, resn='ARG'
alter resi 52 and level20, resn='ALA'
alter resi 87 and level20, resn='THR'
alter resi 102 and level20, resn='SER'
alter resi 107 and level20, resn='SER'
alter resi 142 and level20, resn='SER'
alter resi 181 and level20, resn='SER'
alter resi 80 and level20, resn='TYR'
h_add level20
save mutant_pdbs/GPX6_level20.pdb, level20
delete level20

