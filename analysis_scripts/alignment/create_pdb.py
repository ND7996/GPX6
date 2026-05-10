from Bio.PDB import PDBIO, StructureBuilder
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import numpy as np

sequence = (
    "MIRQLWASSLFPLFLVGFAQLTPESQKMKMDCYKGVTGTIYEYGALTLNGEEYI"
    "QFKQYAGKHVLFINVATYGLTAQYPELNALQEELKPFGVVLLGFPCNQFGKQE"
    "PGKNSEILSGLKYVRPGGGFVPNFQLFEKGDVNGEKEQKVFTFLKNSCPPTSDL"
    "LGSPKQLFWEPMKVHDIRWNFEKFLVGPDGVPVMRWFHRASVSTVKSDILEYLK"
    "QFTPE"
)

aa_map = {
    "A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY",
    "H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","N":"ASN",
    "P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL",
    "W":"TRP","Y":"TYR","U":"SEC"
}

structure = Structure("AF-F7CPY6-F1")
model = Model(0)
chain = Chain("A")

x = 0.0
for i, aa in enumerate(sequence, start=1):
    resname = aa_map[aa]
    residue = Residue((" ", i, " "), resname, "")

    # Idealized backbone geometry (extended chain)
    N  = np.array([x,     0.0, 0.0])
    CA = np.array([x+1.45,0.0, 0.0])
    C  = np.array([x+2.95,0.0, 0.0])
    O  = np.array([x+3.55,0.6, 0.0])

    residue.add(Atom("N",  N,  1.0, 0.0, " ", "N",  None, "N"))
    residue.add(Atom("CA", CA, 1.0, 0.0, " ", "CA", None, "C"))
    residue.add(Atom("C",  C,  1.0, 0.0, " ", "C",  None, "C"))
    residue.add(Atom("O",  O,  1.0, 0.0, " ", "O",  None, "O"))

    chain.add(residue)
    x += 3.8  # Cα–Cα spacing

model.add(chain)
structure.add(model)

io = PDBIO()
io.set_structure(structure)
io.save("AF-F7CPY6-F1_backbone.pdb")

print("Written: AF-F7CPY6-F1_backbone.pdb")
