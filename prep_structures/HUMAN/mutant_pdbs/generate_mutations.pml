# PyMOL script for GPX6 human mutations
# Load this script in PyMOL: File -> Run Script -> select this file

load /home/hp/nayanika/github/GPX6/humanWT/humanWT.pdb, original

# Level 1: T87K
copy level1, original
alter resi 87 and level1, resn='LYS'
h_add level1
save mutant_pdbs/GPX6_level01.pdb, level1
delete level1

# Level 2: T87K, A47S
copy level2, original
alter resi 87 and level2, resn='LYS'
alter resi 47 and level2, resn='SER'
h_add level2
save mutant_pdbs/GPX6_level02.pdb, level2
delete level2

# Level 3: T87K, A47S, S143E
copy level3, original
alter resi 87 and level3, resn='LYS'
alter resi 47 and level3, resn='SER'
alter resi 143 and level3, resn='GLU'
h_add level3
save mutant_pdbs/GPX6_level03.pdb, level3
delete level3

# Level 4: T87K, A47S, S143E, A60T
copy level4, original
alter resi 87 and level4, resn='LYS'
alter resi 47 and level4, resn='SER'
alter resi 143 and level4, resn='GLU'
alter resi 60 and level4, resn='THR'
h_add level4
save mutant_pdbs/GPX6_level04.pdb, level4
delete level4

# Level 5: T87K, A47S, S143E, A60T, F104Y
copy level5, original
alter resi 87 and level5, resn='LYS'
alter resi 47 and level5, resn='SER'
alter resi 143 and level5, resn='GLU'
alter resi 60 and level5, resn='THR'
alter resi 104 and level5, resn='TYR'
h_add level5
save mutant_pdbs/GPX6_level05.pdb, level5
delete level5

# Level 6: T87K, A47S, S143E, A60T, F104Y, S142P
copy level6, original
alter resi 87 and level6, resn='LYS'
alter resi 47 and level6, resn='SER'
alter resi 143 and level6, resn='GLU'
alter resi 60 and level6, resn='THR'
alter resi 104 and level6, resn='TYR'
alter resi 142 and level6, resn='PRO'
h_add level6
save mutant_pdbs/GPX6_level06.pdb, level6
delete level6

# Level 7: T87K, A47S, S143E, A60T, F104Y, S142P, L139F
copy level7, original
alter resi 87 and level7, resn='LYS'
alter resi 47 and level7, resn='SER'
alter resi 143 and level7, resn='GLU'
alter resi 60 and level7, resn='THR'
alter resi 104 and level7, resn='TYR'
alter resi 142 and level7, resn='PRO'
alter resi 139 and level7, resn='PHE'
h_add level7
save mutant_pdbs/GPX6_level07.pdb, level7
delete level7

# Level 8: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R
copy level8, original
alter resi 87 and level8, resn='LYS'
alter resi 47 and level8, resn='SER'
alter resi 143 and level8, resn='GLU'
alter resi 60 and level8, resn='THR'
alter resi 104 and level8, resn='TYR'
alter resi 142 and level8, resn='PRO'
alter resi 139 and level8, resn='PHE'
alter resi 181 and level8, resn='ARG'
h_add level8
save mutant_pdbs/GPX6_level08.pdb, level8
delete level8

# Level 9: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T
copy level9, original
alter resi 87 and level9, resn='LYS'
alter resi 47 and level9, resn='SER'
alter resi 143 and level9, resn='GLU'
alter resi 60 and level9, resn='THR'
alter resi 104 and level9, resn='TYR'
alter resi 142 and level9, resn='PRO'
alter resi 139 and level9, resn='PHE'
alter resi 181 and level9, resn='ARG'
alter resi 52 and level9, resn='THR'
h_add level9
save mutant_pdbs/GPX6_level09.pdb, level9
delete level9

# Level 10: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F
copy level10, original
alter resi 87 and level10, resn='LYS'
alter resi 47 and level10, resn='SER'
alter resi 143 and level10, resn='GLU'
alter resi 60 and level10, resn='THR'
alter resi 104 and level10, resn='TYR'
alter resi 142 and level10, resn='PRO'
alter resi 139 and level10, resn='PHE'
alter resi 181 and level10, resn='ARG'
alter resi 52 and level10, resn='THR'
alter resi 48 and level10, resn='PHE'
h_add level10
save mutant_pdbs/GPX6_level10.pdb, level10
delete level10

# Level 11: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H
copy level11, original
alter resi 87 and level11, resn='LYS'
alter resi 47 and level11, resn='SER'
alter resi 143 and level11, resn='GLU'
alter resi 60 and level11, resn='THR'
alter resi 104 and level11, resn='TYR'
alter resi 142 and level11, resn='PRO'
alter resi 139 and level11, resn='PHE'
alter resi 181 and level11, resn='ARG'
alter resi 52 and level11, resn='THR'
alter resi 48 and level11, resn='PHE'
alter resi 144 and level11, resn='HIS'
h_add level11
save mutant_pdbs/GPX6_level11.pdb, level11
delete level11

# Level 12: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H, Q54T
copy level12, original
alter resi 87 and level12, resn='LYS'
alter resi 47 and level12, resn='SER'
alter resi 143 and level12, resn='GLU'
alter resi 60 and level12, resn='THR'
alter resi 104 and level12, resn='TYR'
alter resi 142 and level12, resn='PRO'
alter resi 139 and level12, resn='PHE'
alter resi 181 and level12, resn='ARG'
alter resi 52 and level12, resn='THR'
alter resi 48 and level12, resn='PHE'
alter resi 144 and level12, resn='HIS'
alter resi 54 and level12, resn='THR'
h_add level12
save mutant_pdbs/GPX6_level12.pdb, level12
delete level12

# Level 13: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H, Q54T, Q177H
copy level13, original
alter resi 87 and level13, resn='LYS'
alter resi 47 and level13, resn='SER'
alter resi 143 and level13, resn='GLU'
alter resi 60 and level13, resn='THR'
alter resi 104 and level13, resn='TYR'
alter resi 142 and level13, resn='PRO'
alter resi 139 and level13, resn='PHE'
alter resi 181 and level13, resn='ARG'
alter resi 52 and level13, resn='THR'
alter resi 48 and level13, resn='PHE'
alter resi 144 and level13, resn='HIS'
alter resi 54 and level13, resn='THR'
alter resi 177 and level13, resn='HIS'
h_add level13
save mutant_pdbs/GPX6_level13.pdb, level13
delete level13

# Level 14: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H, Q54T, Q177H, S102G
copy level14, original
alter resi 87 and level14, resn='LYS'
alter resi 47 and level14, resn='SER'
alter resi 143 and level14, resn='GLU'
alter resi 60 and level14, resn='THR'
alter resi 104 and level14, resn='TYR'
alter resi 142 and level14, resn='PRO'
alter resi 139 and level14, resn='PHE'
alter resi 181 and level14, resn='ARG'
alter resi 52 and level14, resn='THR'
alter resi 48 and level14, resn='PHE'
alter resi 144 and level14, resn='HIS'
alter resi 54 and level14, resn='THR'
alter resi 177 and level14, resn='HIS'
alter resi 102 and level14, resn='GLY'
h_add level14
save mutant_pdbs/GPX6_level14.pdb, level14
delete level14

# Level 15: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H, Q54T, Q177H, S102G, L24I
copy level15, original
alter resi 87 and level15, resn='LYS'
alter resi 47 and level15, resn='SER'
alter resi 143 and level15, resn='GLU'
alter resi 60 and level15, resn='THR'
alter resi 104 and level15, resn='TYR'
alter resi 142 and level15, resn='PRO'
alter resi 139 and level15, resn='PHE'
alter resi 181 and level15, resn='ARG'
alter resi 52 and level15, resn='THR'
alter resi 48 and level15, resn='PHE'
alter resi 144 and level15, resn='HIS'
alter resi 54 and level15, resn='THR'
alter resi 177 and level15, resn='HIS'
alter resi 102 and level15, resn='GLY'
alter resi 24 and level15, resn='ILE'
h_add level15
save mutant_pdbs/GPX6_level15.pdb, level15
delete level15

# Level 16: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H, Q54T, Q177H, S102G, L24I, N3K
copy level16, original
alter resi 87 and level16, resn='LYS'
alter resi 47 and level16, resn='SER'
alter resi 143 and level16, resn='GLU'
alter resi 60 and level16, resn='THR'
alter resi 104 and level16, resn='TYR'
alter resi 142 and level16, resn='PRO'
alter resi 139 and level16, resn='PHE'
alter resi 181 and level16, resn='ARG'
alter resi 52 and level16, resn='THR'
alter resi 48 and level16, resn='PHE'
alter resi 144 and level16, resn='HIS'
alter resi 54 and level16, resn='THR'
alter resi 177 and level16, resn='HIS'
alter resi 102 and level16, resn='GLY'
alter resi 24 and level16, resn='ILE'
alter resi 3 and level16, resn='LYS'
h_add level16
save mutant_pdbs/GPX6_level16.pdb, level16
delete level16

# Level 17: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H, Q54T, Q177H, S102G, L24I, N3K, H173R
copy level17, original
alter resi 87 and level17, resn='LYS'
alter resi 47 and level17, resn='SER'
alter resi 143 and level17, resn='GLU'
alter resi 60 and level17, resn='THR'
alter resi 104 and level17, resn='TYR'
alter resi 142 and level17, resn='PRO'
alter resi 139 and level17, resn='PHE'
alter resi 181 and level17, resn='ARG'
alter resi 52 and level17, resn='THR'
alter resi 48 and level17, resn='PHE'
alter resi 144 and level17, resn='HIS'
alter resi 54 and level17, resn='THR'
alter resi 177 and level17, resn='HIS'
alter resi 102 and level17, resn='GLY'
alter resi 24 and level17, resn='ILE'
alter resi 3 and level17, resn='LYS'
alter resi 173 and level17, resn='ARG'
h_add level17
save mutant_pdbs/GPX6_level17.pdb, level17
delete level17

# Level 18: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H, Q54T, Q177H, S102G, L24I, N3K, H173R, A178T
copy level18, original
alter resi 87 and level18, resn='LYS'
alter resi 47 and level18, resn='SER'
alter resi 143 and level18, resn='GLU'
alter resi 60 and level18, resn='THR'
alter resi 104 and level18, resn='TYR'
alter resi 142 and level18, resn='PRO'
alter resi 139 and level18, resn='PHE'
alter resi 181 and level18, resn='ARG'
alter resi 52 and level18, resn='THR'
alter resi 48 and level18, resn='PHE'
alter resi 144 and level18, resn='HIS'
alter resi 54 and level18, resn='THR'
alter resi 177 and level18, resn='HIS'
alter resi 102 and level18, resn='GLY'
alter resi 24 and level18, resn='ILE'
alter resi 3 and level18, resn='LYS'
alter resi 173 and level18, resn='ARG'
alter resi 178 and level18, resn='THR'
h_add level18
save mutant_pdbs/GPX6_level18.pdb, level18
delete level18

# Level 19: T87K, A47S, S143E, A60T, F104Y, S142P, L139F, S181R, A52T, Y48F, Q144H, Q54T, Q177H, S102G, L24I, N3K, H173R, A178T, A74G
copy level19, original
alter resi 87 and level19, resn='LYS'
alter resi 47 and level19, resn='SER'
alter resi 143 and level19, resn='GLU'
alter resi 60 and level19, resn='THR'
alter resi 104 and level19, resn='TYR'
alter resi 142 and level19, resn='PRO'
alter resi 139 and level19, resn='PHE'
alter resi 181 and level19, resn='ARG'
alter resi 52 and level19, resn='THR'
alter resi 48 and level19, resn='PHE'
alter resi 144 and level19, resn='HIS'
alter resi 54 and level19, resn='THR'
alter resi 177 and level19, resn='HIS'
alter resi 102 and level19, resn='GLY'
alter resi 24 and level19, resn='ILE'
alter resi 3 and level19, resn='LYS'
alter resi 173 and level19, resn='ARG'
alter resi 178 and level19, resn='THR'
alter resi 74 and level19, resn='GLY'
h_add level19
save mutant_pdbs/GPX6_level19.pdb, level19
delete level19

