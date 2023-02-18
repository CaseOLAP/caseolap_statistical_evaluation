import os
import pandas as pd
from IPython.display import display
import xml.etree.ElementTree as ET
#from grape import Graph


def get_answer(msg):
    '''
    FUNCTION:
    - Get a y/n answer from a prompt
    PARAMS:
    - msg (str): The message/question to prompt the user.
    OUTPUT:
    - ans (str): the variable for the answer
    '''
    ans = -1
    while ans != 'y':
        ans = input(msg)
        if ans == 'n':
            break
    return ans
    

def prompt_user_mesh():
    '''
    FUNCTION:
    - Ask the user if they want to create the pre-graph MeSH files
    '''
    # MeSH ID to Tree
    mesh_id_ans = get_answer('Prep the MeSH ID to Tree files? y/n: ')
    
    # MeSH Tree to Tree
    mesh_tree_ans = get_answer('Prep the MeSH Tree to Tree files? y/n: ')

    return mesh_id_ans, mesh_tree_ans
            


def prep_mesh_download(year, input_dir):  
    '''
    FUNCTION:
    - Download the MeSH xml file with IDs and tree numbers
    PARAMS:
    - year (str): The year of the MeSH file xml to download
    - input_dir (str): The directory where input data will be stored
    '''
    year = '2023'
    url = f'https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc{year}.xml'
    os.system(f'wget -N -P {input_dir} {url}')
    if not os.path.exists(input_dir):
        os.system(f'mkdir {input_dir}')
    tree = ET.parse(f'data/desc{year}.xml')
    root = tree.getroot()
    return root



def prep_mesh_id_to_tree(root, directional, output_dir):
    '''
    FUNCTION:
    - This prepares node and edge tables mapping the MeSH Tree Numbers to the 
      MeSH IDs. The MeSH Tree Numbers are part of a tree.
    PARAMS:
    - root (xml-like object): The root of the tree from the MeSH XML
    - directional (bool): indicating if directional edges are used in the graph
    - output_dir (str): The output directory for the MeSH ID to tree mapping
    '''

    '''Get initial data: ID to tree '''
    tree_to_id, id2tree = dict(), dict()
    heads, relations, tails, weights = [],[],[],[]
    rel = '-is-'
    all_tree_numbers = list()
    disease_prefix = ('C','F03')

    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')

            # If disease
            for tree_number in tree_numbers:
                if tree_number.text.startswith(disease_prefix):# and not tree_number.text.startswith('C22'): # exclude animal diseases C22?
                    tree_number = tree_number.text                

                    try:
                        # MeSH ID
                        ID = ele.find('DescriptorUI').text

                        # MeSH ID -[is]- MeSH Tree
                        id2tree.setdefault(ID,set()).add(tree_number)
                        tree_to_id[tree_number] = ID
                        heads.append(tree_number)
                        relations.append(rel)
                        tails.append(ID)
                        weights.append(1.0)
                    except:
                        pass
        except:
            pass

    '''Nodes'''
    tree_numbers = list(tree_to_id.keys())
    mesh_ids = list(id2tree.keys())

    nodes = tree_numbers + mesh_ids
    node_types = ['MeSH_Tree_Disease']*len(tree_numbers) +\
                 ['MeSH_ID_Disease']*len(mesh_ids)
    mesh_tree_to_id_nodes = pd.DataFrame({
                            'node': nodes,
                            'node_type': node_types})

    mesh_tree_to_id_nodes.to_csv(output_dir+'/nodes_mesh_tree_to_id.csv')

    '''Edges'''
    if directional:
        heads += tails
        relations += relations
        tails += heads
        weights += weights
        
    edges = {'head':heads,'relation':relations,'tail':tails,'weight':weights}
    mesh_tree_to_id_edges = pd.DataFrame(edges)
    mesh_tree_to_id_edges.to_csv(output_dir+'/edges_mesh_tree_to_id.csv',
                                index = False)

    
    
def prep_mesh_tree_to_tree(root, output_dir):
    '''
    FUNCTION: 
    - Prepare the MeSH tree data: tree numbers in a hierarchy
    PARAMS:
    - root: tree object from parsed MeSH xml
    - output_dir (str): The output directory for the MeSH tree to tree mapping
    '''
    
    '''Nodes'''
    # Get tree numbers
    all_tree_numbers = list()
    disease_prefix = ('C','F03')
    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')

            # If disease
            for tree_number in tree_numbers:
                if tree_number.text.startswith(disease_prefix):# and not tree_number.text.startswith('C22'): # exclude animal diseases C22?
                    tree_number = tree_number.text
                    all_tree_numbers.append(tree_number)
        except:
            pass
    all_tree_numbers.append('MeSH_Tree_Disease_Root')

    nodes = all_tree_numbers
    node_types = ['MeSH_Tree_Disease']*len(nodes)
    node_dict = {'node': nodes, 'node_type': node_types}
    dis_tree_node_df = pd.DataFrame(node_dict)
    dis_tree_node_df.to_csv(output_dir+'/nodes_mesh_tree_to_tree.csv',
                            index=False)
    
    ''' Edges '''
    tree_to_tree = dict()
    heads, relations, tails, weights = [], [], [], []
    rel = '-MeSH_Tree_Parent_of->'

    # Child
    for tree_num in all_tree_numbers:
        if '.' in tree_num:

            # Parent
            parent = ''
            for num in tree_num.split('.')[:-1]:
                parent += num+'.'
            parent = parent.strip('.')

            # Edges
            heads.append(parent)
            relations.append(rel)
            tails.append(tree_num)
            weights.append(1.0)

    # Add top level to connect trees at the top
    level_one_nodes = {node for node in heads if '.' not in node}
    for level_one_node in level_one_nodes:
        heads.append('MeSH_Tree_Disease_Root')
        relations.append(rel)
        tails.append(level_one_node)
        weights.append(1.0)

    # Create edge df
    dis_tree_edges = pd.DataFrame({'head':heads,
                       'relation':relations,
                       'tail':tails,
                       'weight':weights})
    dis_tree_edges.to_csv(output_dir+'/edges_mesh_tree_to_tree.csv',
                          index=False)


    
def prep_all_mesh_files(mesh_id_ans, mesh_tree_ans, 
                        directional, year='2023', 
                        input_dir = 'data', output_dir='output'):
    '''
    FUNCTION:
    - Prep the MeSH files (tree to tree, tree to ID)
    PARAMS:
    - mesh_id_ans (str): y/n indicating if the MeSH ID to tree file should be created 
    - mesh_tree_ans (str): y/n indicating if the MeSH tree to tree file should be created
    - directional (bool): indicating if directional edges are used in the graph
    - year (str): indicates the year of the MeSH file to download
    - output_dir (str): the output directory where the node and edge files will be stored
    '''
    # Check whether the output directory exists
    if not os.path.exists(output_dir):
        os.system(f'mkdir {output_dir}')
    
    # Download xml?
    if mesh_id_ans == 'y' or mesh_tree_ans == 'y':
        root = prep_mesh_download(year, input_dir)
        
    # Prep ID to tree
    if mesh_id_ans == 'y':
        prep_mesh_id_to_tree(root, directional, output_dir)
        print('Created MeSH ID -is- MeSH Tree files')
    
    # Prep tree to tree
    if mesh_tree_ans == 'y':
        prep_mesh_tree_to_tree(root, output_dir)
        print('Created MeSH Tree -is- MeSH Tree files')
    
    
def merge_all_nodes_and_edges(output_dir):
    '''
    FUNCTION:
    - Take the different node and edge files as input (e.g., MeSH ID to tree), 
      merge and output them into two main csv files: nodes and edges
    '''
    
    '''Nodes'''
    # Import
    node_cols = ['node','node_type']
    nodes_tree_to_id_path = os.path.join(output_dir, 'nodes_mesh_tree_to_tree.csv')
    nodes_tree_to_tree_path = os.path.join(output_dir, 'nodes_mesh_tree_to_id.csv')
    nodes_mesh_tree_to_id = pd.read_csv(nodes_tree_to_id_path)[node_cols]    
    nodes_mesh_tree_to_tree = pd.read_csv(nodes_tree_to_tree_path)[node_cols]
    
    # Export
    node_dfs = [nodes_mesh_tree_to_id, nodes_mesh_tree_to_tree]
    nodes_df = pd.concat(node_dfs).drop_duplicates()
    nodes_df.to_csv(os.path.join(output_dir,'nodes_list.csv'), index=False)

    # Display    
    print('Created nodes dataframe\n')
    display_node_df = get_answer('Do you want to see some nodes? y/n: ')
    if display_node_df == 'y':
        display(nodes_df)
    
    
    '''Edges'''
    weighted_ans = get_answer('Do you want the graph to be weighted? y/n: ')
    
    # Import
    edge_cols = ['head','relation','tail']
    if weighted == True:
        edge_cols.append('weight')
    edges_tree_to_id_path = os.path.join(output_dir,'edges_mesh_tree_to_tree.csv')
    edges_tree_to_tree_path = os.path.join(output_dir,'edges_mesh_tree_to_id.csv')
    edges_mesh_tree_to_id = pd.read_csv(edges_tree_to_id_path)[edge_cols]
    edges_mesh_tree_to_tree = pd.read_csv(edges_tree_to_tree_path)[edge_cols]
    
    # Export
    edges_dfs = [edges_mesh_tree_to_id, edges_mesh_tree_to_tree]
    edges_df = pd.concat(edges_dfs).drop_duplicates()
    edges_df.to_csv(os.path.join(output_dir,'edges_list.csv'), index=False)

    # Display
    print('Created edges dataframe\n')
    display_edges_df = get_answer('Do you want to see some edges? y/n: ')
    if display_edges_df == 'y':
        display(edges_df)

'''Main'''
# Prompt user about prepping the MeSH ID-Tree and Tree-Tree files
mesh_id_ans, mesh_tree_ans = prompt_user_mesh()

# Prep MeSH files
prep_all_mesh_files(mesh_id_ans, mesh_tree_ans, directional=False)
        
# Merge all node files, merge all edge files
merge_all_nodes_and_edges(output_dir = 'output')