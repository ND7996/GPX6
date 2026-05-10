#!/usr/bin/env python3
"""
CREATE LOGICALLY CONNECTED TREE - FIXED VERSION
Connects GPX6 sequences to their corresponding species
"""

import re

def extract_gpx6_subtrees(gpx6_tree):
    """Extract Human_WT and Mouse_WT subtrees from GPX6 tree."""
    print("\n🔍 EXTRACTING GPX6 SUBTREES...")
    
    # First, let's understand the structure better
    # GPX6 tree structure: (Ancestor_node25,((Mouse_WT,...)42,(Human_WT,...)62)63)
    
    # Strategy: Find the closing parenthesis for each subtree
    def find_matching_parenthesis(text, start_pos):
        """Find matching closing parenthesis."""
        count = 1
        pos = start_pos
        while count > 0 and pos < len(text):
            pos += 1
            if text[pos] == '(':
                count += 1
            elif text[pos] == ')':
                count -= 1
        return pos
    
    # Find Mouse_WT subtree
    mouse_start = gpx6_tree.find("(Mouse_WT")
    if mouse_start != -1:
        mouse_end = find_matching_parenthesis(gpx6_tree, mouse_start)
        mouse_subtree = gpx6_tree[mouse_start:mouse_end + 1]  # Include the closing ')'
        print(f"✅ Mouse_WT subtree found: {mouse_subtree[:50]}...")
        
        # Find the node number after Mouse_WT subtree
        # Look for )NN after the closing parenthesis
        match = re.search(r'\)(\d+)', gpx6_tree[mouse_end:mouse_end+10])
        if match:
            mouse_node = match.group(1)
            mouse_subtree = mouse_subtree + mouse_node
            print(f"   Node number: {mouse_node}")
    else:
        print("❌ Could not find Mouse_WT in tree")
        mouse_subtree = None
    
    # Find Human_WT subtree
    human_start = gpx6_tree.find("(Human_WT")
    if human_start != -1:
        human_end = find_matching_parenthesis(gpx6_tree, human_start)
        human_subtree = gpx6_tree[human_start:human_end + 1]
        print(f"✅ Human_WT subtree found: {human_subtree[:50]}...")
        
        # Find the node number
        match = re.search(r'\)(\d+)', gpx6_tree[human_end:human_end+10])
        if match:
            human_node = match.group(1)
            human_subtree = human_subtree + human_node
            print(f"   Node number: {human_node}")
    else:
        print("❌ Could not find Human_WT in tree")
        human_subtree = None
    
    return mouse_subtree, human_subtree

def renumber_subtree(subtree, offset):
    """Renumber nodes in a subtree."""
    if not subtree:
        return None
    
    # Find all node numbers
    nodes = []
    for match in re.finditer(r'\)(\d+)', subtree):
        nodes.append(int(match.group(1)))
    
    # Sort in reverse to avoid replacement conflicts
    nodes = sorted(set(nodes), reverse=True)
    
    # Create a copy to modify
    renumbered = subtree
    
    # First pass: mark nodes for replacement
    for node in nodes:
        renumbered = renumbered.replace(f'){node}', f')TEMP{node}')
    
    # Second pass: replace with new numbers
    for node in nodes:
        renumbered = renumbered.replace(f'TEMP{node}', str(node + offset))
    
    return renumbered

def create_logically_connected_tree():
    """
    Create tree where:
    - Mouse_WT connects to Mus_musculus branch
    - Human_WT connects to Homo_sapiens branch
    """
    
    # Original trees
    gpx6_tree = "(Ancestor_node25,((Mouse_WT,((M2H_Level09,(((M2H_Level10,M2H_Level11)23,(M2H_Level12,M2H_Level13)24)25,((M2H_Level20,(M2H_Level18,M2H_Level19)26)27,((M2H_Level14,M2H_Level15)28,(M2H_Level16,M2H_Level17)29)30)31)32)33,(M2H_Level08,((M2H_Level03,(M2H_Level01,M2H_Level02)34)35,((M2H_Level04,M2H_Level05)36,(M2H_Level06,M2H_Level07)37)38)39)40)41)42,(Human_WT,((H2M_Level19,(H2M_Level17,H2M_Level18)43)44,((((H2M_Level01,H2M_Level02)45,(H2M_Level03,H2M_Level04)46)47,((H2M_Level05,H2M_Level06)48,(H2M_Level07,H2M_Level08)49)50)51,(((H2M_Level09,H2M_Level10)52,(H2M_Level11,H2M_Level12)53)54,((H2M_Level13,H2M_Level14)55,(H2M_Level15,H2M_Level16)56)57)58)59)60)61)62)63"
    
    mammal_tree = "((((((((9_Pan_troglodytes,1_Homo_sapiens)30,12_Macaca_mulatta)29,(3_Callithrix_jacchus,4_Saimiri_boliviensis)31)28,8_Tarsius_syrichta)27,21_Otolemur_garnettii)26,((((((2_Mus_musculus,18_Rattus_norvegicus)37,(11_Mesocricetus_auratus,16_Cricetulus_griseus)38)36,13_Jaculus_jaculus)35,22_Dipodomys_ordii)34,(14_Ictidomys_tridecemlineatus,5_Cavia_porcellus)39)33,10_Oryctolagus_cuniculus)32)25,(((17_Bos_taurus,6_Sus_scrofa)42,7_Equus_caballus)41,(15_Felis_catus,20_Odobenus_rosmarus)43)40)24,19_Loxodonta_africana)23"
    
    print("\n" + "="*80)
    print("CREATING LOGICALLY CONNECTED TREE")
    print("="*80)
    
    # Extract GPX6 subtrees
    mouse_subtree, human_subtree = extract_gpx6_subtrees(gpx6_tree)
    
    if not mouse_subtree or not human_subtree:
        print("\n❌ ERROR: Could not extract GPX6 subtrees")
        print("Using alternative approach...")
        
        # Alternative: Use the entire GPX6 tree but connect properly
        print("\nUsing entire GPX6 tree as reference...")
        
        # Let's manually construct what we need
        mouse_subtree = "(Mouse_WT,((M2H_Level09,(((M2H_Level10,M2H_Level11)23,(M2H_Level12,M2H_Level13)24)25,((M2H_Level20,(M2H_Level18,M2H_Level19)26)27,((M2H_Level14,M2H_Level15)28,(M2H_Level16,M2H_Level17)29)30)31)32)33,(M2H_Level08,((M2H_Level03,(M2H_Level01,M2H_Level02)34)35,((M2H_Level04,M2H_Level05)36,(M2H_Level06,M2H_Level07)37)38)39)40)41)42"
        
        human_subtree = "(Human_WT,((H2M_Level19,(H2M_Level17,H2M_Level18)43)44,((((H2M_Level01,H2M_Level02)45,(H2M_Level03,H2M_Level04)46)47,((H2M_Level05,H2M_Level06)48,(H2M_Level07,H2M_Level08)49)50)51,(((H2M_Level09,H2M_Level10)52,(H2M_Level11,H2M_Level12)53)54,((H2M_Level13,H2M_Level14)55,(H2M_Level15,H2M_Level16)56)57)58)59)60)61)62"
    
    print("\n🔍 ANALYZING MAMMALIAN TREE STRUCTURE:")
    print("  • Homo_sapiens at node 30 with Pan_troglodytes")
    print("  • Mus_musculus at node 37 with Rattus_norvegicus")
    
    # Renumber GPX6 subtrees to avoid conflicts
    offset = 100
    print(f"\n🔄 Renumbering GPX6 nodes with offset: +{offset}")
    
    mouse_gpx6_renum = renumber_subtree(mouse_subtree, offset)
    human_gpx6_renum = renumber_subtree(human_subtree, offset)
    
    print(f"✅ Mouse GPX6 renumbered: {mouse_gpx6_renum[:60]}...")
    print(f"✅ Human GPX6 renumbered: {human_gpx6_renum[:60]}...")
    
    # Insert into mammal tree
    print("\n🔧 INSERTING GPX6 INTO MAMMALIAN TREE...")
    
    # For Human: Replace (9_Pan_troglodytes,1_Homo_sapiens)30
    # with (9_Pan_troglodytes,(1_Homo_sapiens,human_gpx6_renum)200)30
    new_human_node = 200
    human_insertion = f"(9_Pan_troglodytes,(1_Homo_sapiens,{human_gpx6_renum}){new_human_node})30"
    
    mammal_tree_modified = mammal_tree.replace("(9_Pan_troglodytes,1_Homo_sapiens)30", human_insertion)
    print("✅ Inserted Human GPX6 clade")
    
    # For Mouse: Replace (2_Mus_musculus,18_Rattus_norvegicus)37
    # with ((2_Mus_musculus,mouse_gpx6_renum)201,18_Rattus_norvegicus)37
    new_mouse_node = 201
    mouse_insertion = f"((2_Mus_musculus,{mouse_gpx6_renum}){new_mouse_node},18_Rattus_norvegicus)37"
    
    mammal_tree_modified = mammal_tree_modified.replace("(2_Mus_musculus,18_Rattus_norvegicus)37", mouse_insertion)
    print("✅ Inserted Mouse GPX6 clade")
    
    # Add Ancestor_node25 as outgroup
    print("\n🌳 ADDING ANCESTOR AS OUTGROUP...")
    
    # Extract ancestor (it's the first part before the comma)
    ancestor_match = re.match(r'^\(([^,]+),', gpx6_tree)
    if ancestor_match:
        ancestor = ancestor_match.group(1)
        print(f"✅ Found ancestor: {ancestor}")
        
        # Renumber ancestor if needed
        ancestor_renum = ancestor
        if ancestor.endswith('25'):
            ancestor_renum = ancestor + "_renumbered125"
        
        # Create final tree
        root_node = 202
        final_tree = f"({mammal_tree_modified},{ancestor_renum}){root_node};"
    else:
        # Fallback
        root_node = 202
        final_tree = f"({mammal_tree_modified},Ancestor_node25_renumbered125){root_node};"
    
    # Save the tree
    output_file = "LOGICALLY_CONNECTED_TREE_FIXED.nwk"
    
    with open(output_file, 'w') as f:
        f.write(final_tree + "\n")
    
    print(f"\n✅ Created: {output_file}")
    
    # Create a simpler version for clarity
    simple_tree = """
# Simplified logical structure showing connections:
# Tree connects GPX6 sequences to their corresponding species

(
  # Mammalian tree with inserted GPX6 sequences
  (
    ... [other primates] ...,
    
    # Human lineage with experimental sequences
    (9_Pan_troglodytes,
      (1_Homo_sapiens,
        # Human GPX6 experimental evolution
        (Human_WT,
          (H2M_Level19, ... all H2M mutations ...)
        )
      )200
    )30,
    
    ... [more mammals] ...,
    
    # Mouse lineage with experimental sequences  
    (
      (2_Mus_musculus,
        # Mouse GPX6 experimental evolution
        (Mouse_WT,
          (M2H_Level09, ... all M2H mutations ...)
        )
      )201,
      18_Rattus_norvegicus
    )37,
    
    ... [rest of mammals] ...
  )23,
  
  # Ancestral sequence (outgroup)
  Ancestor_node25_renumbered125
)202;
"""
    
    with open("SIMPLIFIED_LOGICAL_TREE.nwk", 'w') as f:
        f.write(simple_tree)
    
    print("✅ Created: SIMPLIFIED_LOGICAL_TREE.nwk")
    
    # Create explanation file
    explanation = """
================================================================
LOGICALLY CONNECTED TREE - EXPLANATION
================================================================

PROBLEM:
In your original merged tree, GPX6 sequences were placed in a 
SEPARATE BRANCH from all mammals. This didn't reflect the 
logical relationship where:
- Human_WT should be close to Homo_sapiens
- Mouse_WT should be close to Mus_musculus

SOLUTION:
We inserted the GPX6 sequences INTO the mammalian tree at the 
appropriate branches:

1. HUMAN LINEAGE:
   Original: (Pan_troglodytes, Homo_sapiens)30
   Modified: (Pan_troglodytes, (Homo_sapiens, Human_GPX6_clade)200)30
   
   Human_GPX6_clade contains:
   - Human_WT (lab version of human GPX6)
   - All H2M mutation levels (evolved toward mouse-like)

2. MOUSE LINEAGE:
   Original: (Mus_musculus, Rattus_norvegicus)37
   Modified: ((Mus_musculus, Mouse_GPX6_clade)201, Rattus_norvegicus)37
   
   Mouse_GPX6_clade contains:
   - Mouse_WT (lab version of mouse GPX6)
   - All M2H mutation levels (evolved toward human-like)

3. ANCESTOR:
   Ancestor_node25 is placed as an outgroup to root the tree.

WHY THIS IS BETTER:
- Phylogenetic relationships match logical relationships
- Mouse_WT is now sister to Mus_musculus (as it should be)
- Human_WT is now sister to Homo_sapiens (as it should be)
- The tree structure reflects the TRUE biological relationships

HOW TO USE:
1. Open LOGICALLY_CONNECTED_TREE_FIXED.nwk in FigTree
2. Compare with your original tree
3. Note how GPX6 sequences are now nested within mammalian branches

================================================================
Node Number Summary:
- Mammalian tree: nodes 23-43
- Human GPX6 clade: nodes 100-162 (original + 100 offset)
- Mouse GPX6 clade: nodes 100-142 (original + 100 offset)
- New connection nodes: 200 (human), 201 (mouse)
- Root: node 202
- Ancestor: renumbered to avoid conflicts
================================================================
"""
    
    with open("TREE_EXPLANATION_FIXED.txt", 'w', encoding='utf-8') as f:
        f.write(explanation)
    
    print("✅ Created: TREE_EXPLANATION_FIXED.txt")
    
    return final_tree

def main():
    """Main function."""
    
    print("\n" + "="*80)
    print("FIXING TREE STRUCTURE - GPX6 Connected to Correct Species")
    print("="*80)
    
    print("""
    This script fixes the problem where GPX6 sequences
    were placed in a separate branch from mammals.
    
    NEW STRUCTURE:
    1. Human_WT connects as sister to Homo_sapiens
    2. Mouse_WT connects as sister to Mus_musculus
    3. This reflects the TRUE logical relationships
    """)
    
    try:
        final_tree = create_logically_connected_tree()
        
        print("\n" + "="*80)
        print("✅ SUCCESS! Tree structure fixed.")
        print("="*80)
        
        print(f"""
        📁 GENERATED FILES:
        
        1. LOGICALLY_CONNECTED_TREE_FIXED.nwk
           - Main tree file for phylogenetic software
           - GPX6 sequences now correctly connected
        
        2. SIMPLIFIED_LOGICAL_TREE.nwk  
           - Simplified version showing the structure
        
        3. TREE_EXPLANATION_FIXED.txt
           - Detailed explanation of changes
        
        🔍 WHAT WAS CHANGED:
        
        BEFORE (Wrong):              AFTER (Correct):
        Root                         Root
        ├── All mammals              ├── Mammals WITH GPX6 inserted
        └── All GPX6                │   ├── (Homo_sapiens + Human_GPX6)
             ├── Human_WT           │   ├── (Mus_musculus + Mouse_GPX6)
             └── Mouse_WT           │   └── ... other mammals
                                   └── Ancestor (outgroup)
        
        🎯 KEY IMPROVEMENT:
        - Human_WT is now PHYLOGENETICALLY CLOSE to Homo_sapiens
        - Mouse_WT is now PHYLOGENETICALLY CLOSE to Mus_musculus
        - Tree structure matches biological reality
        
        💡 NEXT STEPS:
        1. Open LOGICALLY_CONNECTED_TREE_FIXED.nwk in FigTree
        2. Color branches to distinguish:
           - Natural sequences (blue)
           - Experimental sequences (red)
        3. Verify the connections make biological sense
        """)
        
        # Show preview
        print("\n🌳 TREE PREVIEW (first 200 chars):")
        print("  " + final_tree[:200] + "...")
        
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        print("\nTrying alternative approach...")
        
        # Create a minimal working tree
        minimal_tree = "((((((((9_Pan_troglodytes,(1_Homo_sapiens,Human_WT)200)30,12_Macaca_mulatta)29,(3_Callithrix_jacchus,4_Saimiri_boliviensis)31)28,8_Tarsius_syrichta)27,21_Otolemur_garnettii)26,(((((((2_Mus_musculus,Mouse_WT)201,18_Rattus_norvegicus)37,(11_Mesocricetus_auratus,16_Cricetulus_griseus)38)36,13_Jaculus_jaculus)35,22_Dipodomys_ordii)34,(14_Ictidomys_tridecemlineatus,5_Cavia_porcellus)39)33,10_Oryctolagus_cuniculus)32)25,(((17_Bos_taurus,6_Sus_scrofa)42,7_Equus_caballus)41,(15_Felis_catus,20_Odobenus_rosmarus)43)40)24,19_Loxodonta_africana)23,Ancestor_node25)202;"
        
        with open("MINIMAL_LOGICAL_TREE.nwk", 'w') as f:
            f.write(minimal_tree)
        
        print("✅ Created MINIMAL_LOGICAL_TREE.nwk as fallback")
        print("This tree has Human_WT and Mouse_WT connected to their species,")
        print("but without all the mutation levels for simplicity.")

if __name__ == "__main__":
    main()